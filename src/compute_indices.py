# ============================================================
# compute_indices.py
# Compute SPI, SPEI (and the Thornthwaite PET needed for SPEI)
# from HUC4-aggregated monthly time series.
#
# Package: climate_indices (pip install climate_indices)
# ============================================================

import math
import calendar
import numpy as np
import pandas as pd
from climate_indices import indices, compute
from climate_indices.indices import Distribution
from climate_indices.compute import Periodicity


# ============================================================
# Thornthwaite PET
# ============================================================

def _daylight_hours(lat_deg, month):
    """
    Maximum possible daylight hours for a latitude and calendar month.

    Uses the standard astronomical formula (FAO 56, eq. 34):
      ωₛ = arccos(-tan(φ) · tan(δ))   [hour angle at sunset, radians]
      N  = 24/π · ωₛ                   [hours]

    Solar declination δ for the mid-point of each month:
      δ = 0.4093 · sin(2π·J/365 − 1.405)
    where J is the Julian day of the 15th of the month.
    """
    # Approximate Julian day for the 15th of each month
    mid_day = [15, 46, 74, 105, 135, 166, 196, 227, 258, 288, 319, 349]
    J = mid_day[month - 1]

    phi   = math.radians(lat_deg)
    delta = 0.4093 * math.sin(2 * math.pi * J / 365 - 1.405)

    cos_ws = -math.tan(phi) * math.tan(delta)
    cos_ws = max(-1.0, min(1.0, cos_ws))   # clamp to avoid acos domain error
    omega_s = math.acos(cos_ws)
    return 24.0 / math.pi * omega_s        # hours


def thornthwaite_pet(df_temp_c, huc4_lats):
    """
    Compute monthly PET (mm/month) for each HUC4 using Thornthwaite (1948).

    How Thornthwaite PET is calculated:
    ------------------------------------
    Thornthwaite's method derives PET from monthly mean temperature and
    latitude (as a proxy for day length).  It requires no wind, humidity,
    or radiation data — making it ideal when only temperature is available.

    Step 1 — Annual heat index I:
      I = Σₘ max(0, T̄ₘ/5)^1.514   (sum over 12 calendar months)
      where T̄ₘ = long-term mean temperature for month m [°C].
      I integrates the warmth of the year; it is computed once from
      the full time series and held constant.

    Step 2 — Sensitivity exponent a:
      a = 6.75×10⁻⁷ · I³ − 7.71×10⁻⁵ · I² + 1.792×10⁻² · I + 0.49239
      This cubic polynomial was fitted empirically by Thornthwaite.

    Step 3 — Unadjusted monthly PET (for a standardised 30-day, 12 h/day month):
      PET_u = 0                      if T ≤ 0°C   (no evaporation below freezing)
            = 16·(10·T/I)^a          if 0 < T < 26.5°C
            = −415.85 + 32.24·T − 0.43·T²  if T ≥ 26.5°C
      (The piecewise form avoids overflow in (10T/I)^a at very high T.)

    Step 4 — Day-length correction:
      PET_corrected = PET_u × (N_m / 12) × (D_m / 30)
      where N_m = possible sunshine hours for latitude and month (Step 4a),
            D_m = actual days in that calendar month.
      This adjusts from the Thornthwaite standard (12 h, 30 days) to the
      true photoperiod and month length at the HUC4's centroid latitude.

    Accuracy note:
      Thornthwaite PET underestimates in windy, arid climates and
      overestimates in humid climates.  For correlation analysis the
      systematic bias cancels across HUC4s; only the temporal variability
      matters.  Penman-Monteith (used by gridMET) is more accurate but
      requires radiation and wind data.

    Parameters
    ----------
    df_temp_c  : pd.DataFrame, rows=time (monthly DatetimeIndex),
                 cols=HUC4 IDs, values = monthly mean temperature (°C)
    huc4_lats  : dict or pd.Series  {huc4_id: centroid_latitude_degrees_N}

    Returns
    -------
    pd.DataFrame, same shape as df_temp_c, values = PET (mm/month)
    """
    results = {}

    for huc in df_temp_c.columns:
        lat = float(huc4_lats[huc])
        temp = df_temp_c[huc].values.astype(float)

        # Long-term mean temperature for each calendar month (Jan=0, Dec=11)
        # Used to compute the annual heat index I.
        monthly_mean = np.zeros(12)
        times = df_temp_c.index
        for m in range(1, 13):
            mask = times.month == m
            vals = temp[mask]
            monthly_mean[m - 1] = np.nanmean(vals) if vals.size > 0 else 0.0

        # Step 1: Annual heat index I
        I = sum(max(0.0, t / 5) ** 1.514 for t in monthly_mean)
        if I == 0:
            results[huc] = np.zeros(len(temp))
            continue

        # Step 2: Exponent a
        a = 6.75e-7 * I**3 - 7.71e-5 * I**2 + 1.792e-2 * I + 0.49239

        pet_arr = np.zeros(len(temp))
        for i, (t_val, ts) in enumerate(zip(temp, times)):
            if np.isnan(t_val):
                pet_arr[i] = np.nan
                continue

            # Step 3: Unadjusted PET (mm per standard 30-day, 12 h/day month)
            if t_val <= 0:
                pet_u = 0.0
            elif t_val < 26.5:
                pet_u = 16.0 * (10.0 * t_val / I) ** a
            else:
                # High-temperature linear extrapolation (Thornthwaite & Mather 1957)
                pet_u = -415.85 + 32.24 * t_val - 0.43 * t_val**2

            # Step 4: Day-length and month-length correction
            N_m = _daylight_hours(lat, ts.month)
            D_m = calendar.monthrange(ts.year, ts.month)[1]
            pet_arr[i] = pet_u * (N_m / 12.0) * (D_m / 30.0)

        results[huc] = pet_arr

    return pd.DataFrame(results, index=df_temp_c.index)


# ============================================================
# SPI
# ============================================================

def compute_spi(df_precip, scale, calib_start_year, calib_end_year):
    """
    Compute SPI-{scale} for every HUC4 column in df_precip.

    How SPI is calculated (McKee et al. 1993):
    -----------------------------------------
    1. ACCUMULATION
       Roll a {scale}-month window and sum precipitation:
         P_acc[t] = P[t] + P[t-1] + … + P[t-scale+1]   (mm)
       scale=1 → raw monthly total; scale=3 → 3-month total, etc.

    2. FIT A GAMMA DISTRIBUTION (calibration period only)
       The 2-parameter gamma PDF for precipitation x > 0:
         f(x; α, β) = x^(α-1) · exp(-x/β) / (β^α · Γ(α))
       Parameters α (shape) and β (scale) are estimated by maximum
       likelihood, separately for each calendar month (Jan, …, Dec)
       to remove the seasonal cycle.

       Zero-precipitation months are handled by a mixed distribution:
         H(x) = q + (1 - q) · G(x; α, β)
       where q = empirical probability of zero rain in that month.

    3. PROBIT TRANSFORM (standardise to N(0,1))
       SPI = Φ⁻¹(H(P_acc))
       H maps raw precip → uniform probability [0,1];
       Φ⁻¹ (inverse standard normal CDF) converts to a z-score.

    Interpretation:
      SPI > +2.0   extreme wet        SPI < -2.0   extreme drought
      +1 to +2     moderate wet       -2 to -1     moderate drought
      -1 to +1     near normal

    Parameters
    ----------
    df_precip        : pd.DataFrame (mm/month), rows=time, cols=HUC4 IDs
    scale            : int, accumulation window (1 = SPI-1, 3 = SPI-3)
    calib_start_year : int
    calib_end_year   : int

    Returns
    -------
    pd.DataFrame, same shape, values = SPI
    """
    data_start_year = df_precip.index[0].year
    results = {}

    for huc in df_precip.columns:
        arr = df_precip[huc].values.astype(float)
        spi_arr = indices.spi(
            values                   = arr,
            scale                    = scale,
            distribution             = Distribution.gamma,
            data_start_year          = data_start_year,
            calibration_year_initial = calib_start_year,
            calibration_year_final   = calib_end_year,
            periodicity              = Periodicity.monthly,
        )
        results[huc] = spi_arr

    return pd.DataFrame(results, index=df_precip.index)


# ============================================================
# SPEI
# ============================================================

def compute_spei(df_precip, df_pet, scale, calib_start_year, calib_end_year):
    """
    Compute SPEI-{scale} for every HUC4 column.

    How SPEI is calculated (Vicente-Serrano et al. 2010):
    ------------------------------------------------------
    1. WATER BALANCE (monthly deficit / surplus)
         D[t] = P[t] - PET[t]    (mm/month)
       D > 0 → surplus; D < 0 → deficit.
       Unlike SPI, SPEI captures drought driven by high evaporation
       (hot/windy months) even when rainfall is near-average.

    2. ACCUMULATION over {scale} months
         D_acc[t] = D[t] + D[t-1] + … + D[t-scale+1]

    3. FIT A LOG-LOGISTIC DISTRIBUTION (calibration period)
       3-parameter log-logistic PDF:
         f(x; γ, α, β) = (β/α)·((x-γ)/α)^(β-1) / (1+((x-γ)/α)^β)²
       Parameters estimated via L-moments (probability-weighted moments),
       which are more robust than MLE for small samples.
       The climate_indices package uses Pearson III as a close
       approximation to the log-logistic — both fit D_acc similarly.

    4. PROBIT TRANSFORM
         SPEI = Φ⁻¹(F(D_acc))
       Identical to SPI step 3; interpretation is the same.

    Parameters
    ----------
    df_precip        : pd.DataFrame (mm/month)
    df_pet           : pd.DataFrame (mm/month), same shape
    scale            : int (1 or 3)
    calib_start_year : int
    calib_end_year   : int

    Returns
    -------
    pd.DataFrame, same shape, values = SPEI
    """
    common_huc = df_precip.columns.intersection(df_pet.columns)
    common_idx = df_precip.index.intersection(df_pet.index)
    df_precip  = df_precip.loc[common_idx, common_huc]
    df_pet     = df_pet.loc[common_idx, common_huc]

    data_start_year = common_idx[0].year
    results = {}

    for huc in common_huc:
        p_arr   = df_precip[huc].values.astype(float)
        pet_arr = df_pet[huc].values.astype(float)
        spei_arr = indices.spei(
            precips_mm               = p_arr,
            pet_mm                   = pet_arr,
            scale                    = scale,
            distribution             = Distribution.pearson,   # Pearson III ≈ log-logistic
            data_start_year          = data_start_year,
            calibration_year_initial = calib_start_year,
            calibration_year_final   = calib_end_year,
            periodicity              = Periodicity.monthly,
        )
        results[huc] = spei_arr

    return pd.DataFrame(results, index=common_idx)
