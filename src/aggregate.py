# ============================================================
# aggregate.py
# Fetch monthly climate variables from NOAA PSL via OPeNDAP
# and aggregate spatially to HUC4 watershed boundaries.
#
# Why NOAA PSL instead of gridMET:
#   gridMET stores DAILY data (365 × 585 × 1386 cells/year ≈ 1.2 GB
#   per year).  The THREDDS server actively resets the connection
#   when a single request exceeds its transfer limit.
#   NOAA PSL provides MONTHLY datasets (1 value/month vs 30):
#     - CMAP precipitation : 2.5°, already monthly → ~0.3 MB for US/20 yr
#     - Dai PDSI           : 2.5°, already monthly → ~0.3 MB for US/20 yr
#     - GHCN+CAMS temp     : 0.5°, already monthly → ~1 MB for US/20 yr
#   Monthly data is all we need for SPI / SPEI / PDSI correlation
#   and avoids the daily→monthly resampling step entirely.
#
# Why pydap engine:
#   The netcdf4 engine's OPeNDAP client depends on libcurl, which
#   commonly fails on Windows ("curl error details").  pydap is a
#   pure-Python OPeNDAP client that uses the requests library instead.
# ============================================================

import calendar
import numpy as np
import pandas as pd
import xarray as xr
import geopandas as gpd
import regionmask

# --- NOAA PSL OPeNDAP URLs ---
# All datasets are on NOAA's THREDDS server; requests complete in seconds
# because the US monthly slice is tiny (< 2 MB per variable).
NOAA_PSL = {
    # CMAP Enhanced monthly precipitation (mm/day), 2.5°, 1979–present
    # Reference: Xie & Arkin (1997)
    "precip": "https://psl.noaa.gov/thredds/dodsC/Datasets/cmap/enh/precip.mon.mean.nc",

    # Dai self-calibrated PDSI (dimensionless), 2.5°, 1870–2018
    # Reference: Dai (2011, J. Hydrometeorology)
    "pdsi":   "https://psl.noaa.gov/thredds/dodsC/Datasets/dai_pdsi/pdsi.mon.mean.nc",

    # GHCN + CAMS monthly surface air temperature (°C), 0.5°, 1948–present
    # Used to compute Thornthwaite PET for SPEI.
    "temp":   "https://psl.noaa.gov/thredds/dodsC/Datasets/ghcncams/air.mon.mean.nc",
}

# US CONUS bounding box (applied AFTER coordinate normalisation)
US_LAT = (24.0, 50.0)
US_LON = (-130.0, -60.0)


# ============================================================
# Internal helpers
# ============================================================

def _open_psl(url, var_hint):
    """
    Open a NOAA PSL OPeNDAP dataset with pydap, return the target DataArray.

    pydap is tried first (pure Python, no libcurl); netcdf4 is the fallback.
    Prints actual variable / dimension names so mismatches are visible.
    """
    ds = None
    for engine in ("pydap", "netcdf4"):
        try:
            ds = xr.open_dataset(url, engine=engine)
            print(f"  Connected ({engine}) : {url.split('/')[-1]}")
            break
        except Exception as exc:
            print(f"  engine='{engine}' failed: {exc}")

    if ds is None:
        raise ConnectionError(f"Cannot open {url}. Check internet and NOAA PSL status.")

    print(f"  Variables  : {list(ds.data_vars)}")
    print(f"  Dimensions : {dict(ds.sizes)}")

    match = [v for v in ds.data_vars if var_hint.lower() in v.lower()]
    if not match:
        raise KeyError(
            f"No variable containing '{var_hint}' found. "
            f"Available: {list(ds.data_vars)}"
        )
    print(f"  Using      : '{match[0]}'")
    return ds[match[0]]


def _normalize_coords(da):
    """
    Standardise coordinates so slicing always works the same way:
      1. If lon is in 0–360 convention, convert to -180 – 180.
      2. Sort lat ascending (some datasets store 90→-90).
      3. Sort lon ascending.
    This handles CMAP (lon 0–360, lat N→S) and Dai PDSI in one shot.
    """
    # Longitude: 0-360 → -180 to 180
    if float(da.lon.max()) > 180:
        da = da.assign_coords(lon=((da.lon + 180) % 360) - 180)

    da = da.sortby("lat").sortby("lon")
    return da


def _slice_us(da, start, end):
    """Slice DataArray to US extent and requested time window, then download."""
    da = _normalize_coords(da)
    lat_lo, lat_hi = US_LAT
    lon_lo, lon_hi = US_LON

    da_us = da.sel(
        lat=slice(lat_lo, lat_hi),
        lon=slice(lon_lo, lon_hi),
        time=slice(start, end),
    )
    print(f"  Downloading slice {start}–{end} "
          f"({da_us.sizes.get('time','?')} months × "
          f"{da_us.sizes.get('lat','?')} lat × "
          f"{da_us.sizes.get('lon','?')} lon) ...", end=" ", flush=True)
    da_us = da_us.load()
    print("done")
    return da_us


# ============================================================
# Public fetch functions
# ============================================================

def fetch_monthly_precip(start="2000-01", end="2020-12"):
    """
    Pull CMAP monthly precipitation and convert to mm/month.

    Calculation:
      CMAP stores precipitation as a monthly MEAN RATE (mm/day).
      To get a monthly TOTAL comparable to what SPI/SPEI expect:
        P_total [mm/month] = P_rate [mm/day] × days_in_month

      This accounts for the fact that February (28 days) accumulates
      less water than January (31 days) at the same daily rate.

    Returns
    -------
    xr.DataArray  dims = (time, lat, lon),  units = mm/month
    """
    print("Fetching CMAP monthly precipitation from NOAA PSL ...")
    da = _open_psl(NOAA_PSL["precip"], "precip")
    da_us = _slice_us(da, start, end)

    # Multiply mm/day by days_in_month → mm/month
    # time.dt.days_in_month gives the correct value for every month
    days = da_us.time.dt.days_in_month
    da_mm = da_us * days
    da_mm.attrs["units"] = "mm/month"
    return da_mm


def fetch_pdsi(start="2000-01", end="2020-12"):
    """
    Pull Dai self-calibrated PDSI from NOAA PSL.

    How PDSI is calculated (Palmer 1965, self-calibrated variant):
      Palmer's model tracks soil moisture in two layers (surface + deep)
      using a monthly water balance:
        P - ET_actual - R - RO = ΔSoil moisture
      where R = recharge, RO = runoff.

      A climatically-appropriate precipitation P̂ is estimated from
      the long-term mean of each flux component.  The difference
      (P - P̂) is converted to a standardised moisture anomaly Z:
        Z = k · (P - P̂)
      where k is a climate-weighting factor (re-derived per grid cell
      in the self-calibrated version, making regional PDSI comparable).

      PDSI is a recursive weighted sum of monthly Z scores:
        PDSIₜ = 0.897·PDSIₜ₋₁ + Z/3

      The 0.897 multiplier gives PDSI a memory of ~3–6 months,
      so it responds more slowly to rainfall than SPI-1/SPEI-1.
      Positive = wetter than normal; negative = drought.

    Note: Dai (2011) PDSI covers 1870–2018.  For 2000–2018 we have
    19 years; the PDSI column in the correlation summary will reflect
    this shorter window.

    Returns
    -------
    xr.DataArray  dims = (time, lat, lon),  dimensionless
    """
    print("Fetching Dai scPDSI from NOAA PSL ...")
    da = _open_psl(NOAA_PSL["pdsi"], "pdsi")
    da_us = _slice_us(da, start, end)
    return da_us


def fetch_monthly_temp(start="2000-01", end="2020-12"):
    """
    Pull GHCN+CAMS monthly mean surface air temperature (°C).

    Temperature is used to compute Thornthwaite PET for SPEI:
      PET = f(monthly_mean_T [°C], latitude [°])

    Returns
    -------
    xr.DataArray  dims = (time, lat, lon),  units = °C
    """
    print("Fetching GHCN+CAMS monthly temperature from NOAA PSL ...")
    da = _open_psl(NOAA_PSL["temp"], "air")
    da_us = _slice_us(da, start, end)
    return da_us


# ============================================================
# Spatial aggregation
# ============================================================

def aggregate_to_huc4(da, huc4_gdf, huc4_id_col="huc4"):
    """
    Spatially average a gridded DataArray over each HUC4 polygon.

    How the aggregation works (regionmask method):
      Step 1 — Build a 2-D integer mask (lat × lon):
        For every grid cell, assign the row-index of the HUC4 polygon
        whose interior contains the cell centre point.
        Cells outside all HUC4 polygons get NaN.
        Result: mask[lat, lon] ∈ {0, 1, …, N-1, NaN}

      Step 2 — GroupBy + mean:
        da.groupby(mask).mean() collects all cells with the same mask
        value (= same HUC4) and averages them.
        For dims (time, lat, lon) this gives (time, n_huc4).
        NaN cells are automatically excluded.

      Arithmetic mean of grid cells = area-weighted average when cells
      are equal-area.  At 2.5° latitude, cells near 50°N are ~15% smaller
      than cells near 24°N, so strictly speaking a cosine area-weight
      should be applied.  For HUC4-scale correlation analysis the
      difference is small and we use the simpler unweighted mean.

    Parameters
    ----------
    da          : xr.DataArray  dims = (time, lat, lon)
    huc4_gdf    : GeoDataFrame, CRS = EPSG:4326
    huc4_id_col : column holding the HUC4 ID string (default 'huc4')

    Returns
    -------
    pd.DataFrame  index = pd.DatetimeIndex (monthly),  columns = HUC4 IDs
    """
    print(f"Building HUC4 mask on {da.sizes.get('lat','?')}×"
          f"{da.sizes.get('lon','?')} grid ({len(huc4_gdf)} polygons) ...")

    mask = regionmask.mask_geopandas(huc4_gdf, da.lon.values, da.lat.values)
    mask.name = "region"

    print("Computing spatial mean per HUC4 ...")
    grouped = da.groupby(mask).mean()

    # Drop -1 group (older regionmask used -1 for cells outside all polygons)
    if -1 in grouped.region.values:
        grouped = grouped.sel(region=grouped.region >= 0)

    # Map integer row-indices → HUC4 ID strings
    row_indices = grouped.region.values.astype(int)
    huc4_ids    = huc4_gdf[huc4_id_col].values[row_indices]
    grouped     = grouped.assign_coords(region=huc4_ids).rename({"region": "huc4"})

    df = grouped.to_dataframe(name="value").reset_index()
    df = df.pivot(index="time", columns="huc4", values="value")
    df.index = pd.to_datetime(df.index)
    df.columns.name = None
    return df
