# ============================================================
# drought_index_analysis.py
# Main runner: correlate SPI / SPEI / PDSI with rainfall
# at the HUC4 watershed scale for the CONUS, 2000-2020.
#
# Data sources (all via NOAA PSL OPeNDAP — monthly, tiny requests):
#   Precipitation : CMAP Enhanced, 2.5°, 1979-present
#   PDSI          : Dai scPDSI,    2.5°, 1870-2018
#   Temperature   : GHCN+CAMS,    0.5°, 1948-present
#   PET           : Thornthwaite formula from temperature + latitude
#
# SPI and SPEI are computed from the CMAP precipitation series
# using the climate_indices package.
# ============================================================

import os
import pandas as pd
import geopandas as gpd
import matplotlib.pyplot as plt

from src.aggregate       import (fetch_monthly_precip, fetch_pdsi,
                                  fetch_monthly_temp, aggregate_to_huc4)
from src.compute_indices import compute_spi, compute_spei, thornthwaite_pet
from src.correlation     import compute_correlation
from src.plot            import plot_correlation_map, plot_correlation_boxplot

# ============================================================
# Paths
# ============================================================
BASE_DIR = r"C:\Users\anyiw\OneDrive - University of Illinois - Urbana\GIS\2026REU\WSDI"
HUC4_SHP = os.path.join(BASE_DIR, "data", "HUC4", "WBDHU4.shp")
OUT_DIR  = os.path.join(BASE_DIR, "output", "drought_correlation")
FIG_DIR  = os.path.join(OUT_DIR, "figures")
os.makedirs(OUT_DIR, exist_ok=True)
os.makedirs(FIG_DIR, exist_ok=True)

# Time window: aligned with existing WSDI data (2000-2020).
# Note: Dai PDSI ends 2018, so PDSI correlation uses 2000-2018 (19 yr).
START, END             = "2000-01", "2020-12"
CALIB_START, CALIB_END = 2000, 2020

# ============================================================
# STEP 1: Load HUC4 boundary (CONUS only)
# ============================================================
print("=" * 55)
print("STEP 1: Loading HUC4 shapefile ...")
print("=" * 55)

huc4 = gpd.read_file(HUC4_SHP)
huc4 = huc4.to_crs("EPSG:4326")
huc4["huc4"] = huc4["huc4"].astype(str).str.zfill(4)

# HUC4 codes 01-18 = CONUS; 19 = Alaska, 20 = Hawaii, 21 = Puerto Rico
huc4 = huc4[huc4["huc4"].str[:2].astype(int) <= 18].copy().reset_index(drop=True)
print(f"  CONUS HUC4 units: {len(huc4)}")

# Centroid latitude for each HUC4 (needed for Thornthwaite PET)
huc4_lats = dict(zip(
    huc4["huc4"].values,
    huc4.geometry.centroid.y.values
))

# ============================================================
# STEP 2: Precipitation -> aggregate to HUC4
# ============================================================
print("\n" + "=" * 55)
print("STEP 2: CMAP monthly precipitation ...")
print("=" * 55)

da_pr   = fetch_monthly_precip(START, END)
df_rain = aggregate_to_huc4(da_pr, huc4)
df_rain.to_csv(os.path.join(OUT_DIR, "huc4_rainfall_mm.csv"))
print(f"  Shape: {df_rain.shape}  (months × HUC4s)")

# ============================================================
# STEP 3: Temperature -> aggregate to HUC4 -> Thornthwaite PET
# ============================================================
print("\n" + "=" * 55)
print("STEP 3: GHCN+CAMS temperature -> Thornthwaite PET ...")
print("=" * 55)

da_temp  = fetch_monthly_temp(START, END)
df_temp  = aggregate_to_huc4(da_temp, huc4)
df_temp.to_csv(os.path.join(OUT_DIR, "huc4_temp_K.csv"))

# GHCN+CAMS stores temperature in Kelvin; Thornthwaite requires Celsius.
df_temp_c = df_temp - 273.15
df_temp_c.to_csv(os.path.join(OUT_DIR, "huc4_temp_c.csv"))

# Thornthwaite PET: computed from temperature + HUC4 centroid latitude.
# No additional download needed — pure arithmetic on the HUC4 DataFrame.
print("  Computing Thornthwaite PET from HUC4 temperature ...")
df_pet = thornthwaite_pet(df_temp_c, huc4_lats)
df_pet.to_csv(os.path.join(OUT_DIR, "huc4_pet_mm.csv"))
print(f"  PET shape: {df_pet.shape}")

# ============================================================
# STEP 4: PDSI -> aggregate to HUC4
# ============================================================
print("\n" + "=" * 55)
print("STEP 4: Dai scPDSI ...")
print("=" * 55)

da_pdsi = fetch_pdsi(START, END)   # may end at 2018 if Dai file ends there
df_pdsi = aggregate_to_huc4(da_pdsi, huc4)
df_pdsi.to_csv(os.path.join(OUT_DIR, "huc4_pdsi.csv"))
print(f"  PDSI shape: {df_pdsi.shape}  (note: Dai dataset ends ~2018)")

# ============================================================
# STEP 5: Compute SPI-1, SPI-3, SPEI-1, SPEI-3
# ============================================================
print("\n" + "=" * 55)
print("STEP 5: Computing SPI and SPEI ...")
print("=" * 55)

print("  SPI-1  (gamma, 1-month accumulation) ...")
df_spi1 = compute_spi(df_rain, scale=1,
                      calib_start_year=CALIB_START, calib_end_year=CALIB_END)
df_spi1.to_csv(os.path.join(OUT_DIR, "huc4_spi1.csv"))

print("  SPI-3  (gamma, 3-month accumulation) ...")
df_spi3 = compute_spi(df_rain, scale=3,
                      calib_start_year=CALIB_START, calib_end_year=CALIB_END)
df_spi3.to_csv(os.path.join(OUT_DIR, "huc4_spi3.csv"))

print("  SPEI-1 (Pearson III on D = P - PET, 1-month) ...")
df_spei1 = compute_spei(df_rain, df_pet, scale=1,
                        calib_start_year=CALIB_START, calib_end_year=CALIB_END)
df_spei1.to_csv(os.path.join(OUT_DIR, "huc4_spei1.csv"))

print("  SPEI-3 (3-month water balance accumulation) ...")
df_spei3 = compute_spei(df_rain, df_pet, scale=3,
                        calib_start_year=CALIB_START, calib_end_year=CALIB_END)
df_spei3.to_csv(os.path.join(OUT_DIR, "huc4_spei3.csv"))
print("  All indices computed.")

# ============================================================
# STEP 6: Pearson correlation — each index vs same-month rainfall
#
# Why same-month rainfall:
#   SPI-1 and SPEI-1 are 1-month anomalies -> should correlate strongly
#   with same-month rain.  SPI-3/SPEI-3 accumulate 3 months -> a 3-month
#   rainfall average would give a fairer comparison, but we use same-month
#   for consistency across all indices.
#   PDSI has a 3-6 month lag (recursive 0.897 weight) so its same-month r
#   will be lower than SPI-1/SPEI-1 by design.
# ============================================================
print("\n" + "=" * 55)
print("STEP 6: Pearson correlations (index vs same-month rainfall) ...")
print("=" * 55)

indices_map = {
    "SPI-1":  df_spi1,
    "SPI-3":  df_spi3,
    "SPEI-1": df_spei1,
    "SPEI-3": df_spei3,
    "PDSI":   df_pdsi,
}

corr_results = {}
for name, df_idx in indices_map.items():
    corr = compute_correlation(df_idx, df_rain, method="pearson")
    corr_results[name] = corr
    safe = name.lower().replace("-", "")
    corr.to_csv(os.path.join(OUT_DIR, f"corr_{safe}.csv"))
    sig  = (corr["p"] < 0.05).sum()
    print(f"  {name:8s} | median r = {corr['r'].median():.3f}"
          f" | mean r = {corr['r'].mean():.3f}"
          f" | sig (p<0.05): {sig}/{len(corr)} HUC4s")

# ============================================================
# STEP 7: Figures + summary
# ============================================================
print("\n" + "=" * 55)
print("STEP 7: Saving summary and figures ...")
print("=" * 55)

# r-value summary (one column per index, one row per HUC4)
summary = pd.DataFrame({name: corr["r"] for name, corr in corr_results.items()})
summary.to_csv(os.path.join(OUT_DIR, "correlation_summary.csv"))

# p-value summary (same layout as correlation_summary.csv)
p_summary = pd.DataFrame({name: corr["p"] for name, corr in corr_results.items()})
p_summary.to_csv(os.path.join(OUT_DIR, "correlation_p_summary.csv"))
print("  p-value summary saved.")

# ── Stats table (Table 2 in the paper) ────────────────────────────────────────
import numpy as np
rows = []
for name, corr in corr_results.items():
    r = corr["r"].dropna()
    p = corr["p"].dropna()
    n_total = len(r)
    n_sig   = (p < 0.05).sum()
    rows.append({
        "Index":           name,
        "Median_r":        round(r.median(), 3),
        "Mean_r":          round(r.mean(),   3),
        "Std_r":           round(r.std(),    3),
        "N_HUC4":          n_total,
        "N_Significant":   int(n_sig),
        "Pct_Significant": f"{100 * n_sig / n_total:.1f}%",
    })
stats_table = pd.DataFrame(rows)
stats_table.to_csv(os.path.join(OUT_DIR, "stats_table.csv"), index=False)
print("  Stats table (Table 2):")
print(stats_table.to_string(index=False))

plot_correlation_boxplot(
    summary,
    save_path=os.path.join(FIG_DIR, "boxplot_correlation.png")
)
print("  Boxplot saved.")

fig, axes = plt.subplots(2, 3, figsize=(22, 12))
axes = axes.flatten()
for i, (name, corr) in enumerate(corr_results.items()):
    plot_correlation_map(huc4, corr,
                         title=f"{name} vs Same-Month Rainfall  (Pearson r, 2000-2020)",
                         ax=axes[i])
axes[-1].set_visible(False)
plt.tight_layout()
plt.savefig(os.path.join(FIG_DIR, "map_correlation_all.png"), dpi=150, bbox_inches="tight")
print("  Spatial map (r) saved.")

# ── Significance map: which HUC4s are NOT significant (p ≥ 0.05) ─────────────
fig2, axes2 = plt.subplots(2, 3, figsize=(22, 12))
axes2 = axes2.flatten()
for i, (name, corr) in enumerate(corr_results.items()):
    ax = axes2[i]
    merged = huc4.merge(corr[["p"]].reset_index(), on="huc4", how="left")
    merged["sig"] = merged["p"] < 0.05
    merged[merged["sig"]].plot(ax=ax, color="#2196F3", linewidth=0.3, label="p < 0.05")
    merged[~merged["sig"]].plot(ax=ax, color="#F44336", linewidth=0.3, label="p ≥ 0.05")
    ax.set_title(f"{name} — Significance (α = 0.05)", fontsize=10)
    ax.axis("off")
    from matplotlib.patches import Patch
    ax.legend(handles=[Patch(color="#2196F3", label="p < 0.05 (significant)"),
                        Patch(color="#F44336", label="p ≥ 0.05 (not significant)")],
              loc="lower left", fontsize=7)
axes2[-1].set_visible(False)
plt.suptitle("Pearson r Significance Map — All Drought Indices vs Same-Month Rainfall",
             fontsize=12, fontweight="bold")
plt.tight_layout()
plt.savefig(os.path.join(FIG_DIR, "map_significance_all.png"), dpi=150, bbox_inches="tight")
print("  Significance map saved.")

print("\n" + "=" * 55)
print("DONE  ->", OUT_DIR)
print("=" * 55)
print("\nMedian Pearson r  (higher = better tracks rainfall):")
print(summary.median().sort_values(ascending=False).round(3).to_string())
print("\nExpected ranking: SPI-1 ~= SPEI-1 > SPI-3 ~= SPEI-3 > PDSI")
print("(SPI-1/SPEI-1 computed from same precip series; PDSI has multi-month lag)")
