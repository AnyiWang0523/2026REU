# ============================================================
# Correlation Analysis: Discharge vs Water Supply
# Only HUC4s with large reservoirs (missing_ratio < 0.05)
# Analysis: discharge_mean vs PS/IR SW/GW Supply mean
# ============================================================

import geopandas as gpd
import pandas as pd
import numpy as np
from scipy import stats
import os

# ============================================================
# Path Configuration
# ============================================================

BASE_DIR   = r"C:\Users\anyiw\OneDrive - University of Illinois - Urbana\GIS\2026REU\WSDI"
HUC4_SHP   = os.path.join(BASE_DIR, "data", "HUC4", "WBDHU4.shp")
RESV_CSV   = os.path.join(BASE_DIR, "data", "CSV", "GDROM_Resv.csv")
GAGE_CSV   = os.path.join(BASE_DIR, "data", "CSV", "USGS_Gages.csv")
OUTPUT_DIR = os.path.join(BASE_DIR, "output", "correlation")
os.makedirs(OUTPUT_DIR, exist_ok=True)

# Supply CSVs
PS_GW_CSV = os.path.join(BASE_DIR, "output", "PS_GW", "WSDI_PS_GW_HUC4_2000_2020.csv")
PS_SW_CSV = os.path.join(BASE_DIR, "output", "PS_SW", "WSDI_PS_SW_HUC4_2000_2020.csv")
IR_GW_CSV = os.path.join(BASE_DIR, "output", "IR_GW", "WSDI_IR_GW_HUC4_2000_2020.csv")
IR_SW_CSV = os.path.join(BASE_DIR, "output", "IR_SW", "WSDI_IR_SW_HUC4_2000_2020.csv")

# ============================================================
# STEP 1: Load and filter reservoirs (missing_ratio < 0.05)
# ============================================================
print("=" * 50)
print("STEP 1: Loading reservoir data...")
print("=" * 50)

resv = pd.read_csv(RESV_CSV)
resv["missing_ratio"] = pd.to_numeric(resv["missing_ratio"], errors="coerce")
resv = resv[resv["missing_ratio"] < 0.05].copy()
resv = resv.dropna(subset=["lattitude", "longtitude"])
print(f"✓ Valid reservoirs: {len(resv)}")

# Convert to GeoDataFrame
resv_gdf = gpd.GeoDataFrame(
    resv,
    geometry=gpd.points_from_xy(resv["longtitude"], resv["lattitude"]),
    crs="EPSG:4326"
)

# ============================================================
# STEP 2: Load HUC4 and spatial join reservoirs -> HUC4
# ============================================================
print("\n" + "=" * 50)
print("STEP 2: Spatial join reservoirs to HUC4...")
print("=" * 50)

huc4_gdf = gpd.read_file(HUC4_SHP)
huc4_gdf["huc4"] = huc4_gdf["huc4"].astype(str).str.zfill(4)
huc4_gdf = huc4_gdf.to_crs("EPSG:4326")

# Spatial join: which HUC4 does each reservoir fall in?
resv_joined = gpd.sjoin(
    resv_gdf,
    huc4_gdf[["huc4", "geometry"]],
    how="left",
    predicate="within"
)

# HUC4s that have at least one valid reservoir
huc4_with_resv = resv_joined["huc4"].dropna().unique()
print(f"✓ HUC4s with large reservoirs: {len(huc4_with_resv)}")

# Reservoir count and max storage per HUC4
resv_per_huc4 = (
    resv_joined.groupby("huc4")
    .agg(
        resv_count        = ("GRanD_ID", "count"),
        max_storage_af    = ("max_historical_Storage_af", "max"),
        mean_storage_af   = ("max_historical_Storage_af", "mean")
    )
    .reset_index()
    .rename(columns={"huc4": "HUC4"})
)

# ============================================================
# STEP 3: Load and filter discharge data (discharge > 0)
# ============================================================
print("\n" + "=" * 50)
print("STEP 3: Loading discharge data...")
print("=" * 50)

gage = pd.read_csv(GAGE_CSV)
gage = gage[gage["discharge_mean"] > 0].copy()
gage = gage.dropna(subset=["dec_lat_va", "dec_long_va"])
print(f"✓ Valid gage stations: {len(gage)}")

# Convert to GeoDataFrame
gage_gdf = gpd.GeoDataFrame(
    gage,
    geometry=gpd.points_from_xy(gage["dec_long_va"], gage["dec_lat_va"]),
    crs="EPSG:4326"
)

# Spatial join gages -> HUC4
gage_joined = gpd.sjoin(
    gage_gdf,
    huc4_gdf[["huc4", "geometry"]],
    how="left",
    predicate="within"
)

# Keep only HUC4s that have reservoirs
gage_joined = gage_joined[gage_joined["huc4"].isin(huc4_with_resv)]

# Mean discharge per HUC4
discharge_huc4 = (
    gage_joined.groupby("huc4")["discharge_mean"]
    .mean()
    .reset_index()
    .rename(columns={"huc4": "HUC4", "discharge_mean": "Discharge_mean"})
)

print(f"✓ HUC4s with discharge data: {len(discharge_huc4)}")

# ============================================================
# STEP 4: Load Supply data and compute mean (2000-2020)
# ============================================================
print("\n" + "=" * 50)
print("STEP 4: Loading supply data...")
print("=" * 50)

def load_supply_mean(csv_path, supply_col, label):
    df = pd.read_csv(csv_path, dtype={"HUC4": str})
    df["HUC4"] = df["HUC4"].str.zfill(4)
    # Keep only HUC4s with reservoirs
    df = df[df["HUC4"].isin(huc4_with_resv)]
    mean_df = (
        df.groupby("HUC4")[supply_col]
        .mean()
        .reset_index()
        .rename(columns={supply_col: label})
    )
    return mean_df

ps_gw = load_supply_mean(PS_GW_CSV, "GW_Supply", "PS_GW_mean")
ps_sw = load_supply_mean(PS_SW_CSV, "SW_Supply", "PS_SW_mean")
ir_gw = load_supply_mean(IR_GW_CSV, "GW_Supply", "IR_GW_mean")
ir_sw = load_supply_mean(IR_SW_CSV, "SW_Supply", "IR_SW_mean")

print(f"✓ PS GW: {len(ps_gw)} HUC4s")
print(f"✓ PS SW: {len(ps_sw)} HUC4s")
print(f"✓ IR GW: {len(ir_gw)} HUC4s")
print(f"✓ IR SW: {len(ir_sw)} HUC4s")

# ============================================================
# STEP 5: Merge all data
# ============================================================
print("\n" + "=" * 50)
print("STEP 5: Merging all data...")
print("=" * 50)

merged = resv_per_huc4.copy()
for df in [discharge_huc4, ps_gw, ps_sw, ir_gw, ir_sw]:
    merged = merged.merge(df, on="HUC4", how="left")

merged = merged.dropna(subset=["Discharge_mean"])
print(f"✓ Final dataset: {len(merged)} HUC4s")
print(merged.head(5).to_string(index=False))

# ============================================================
# STEP 6: Correlation Analysis
# ============================================================
print("\n" + "=" * 50)
print("STEP 6: Pearson Correlation Analysis...")
print("=" * 50)

supply_cols = ["PS_GW_mean", "PS_SW_mean", "IR_GW_mean", "IR_SW_mean"]
corr_results = []

for col in supply_cols:
    sub = merged[["Discharge_mean", col]].dropna()
    if len(sub) < 5:
        continue
    r, p = stats.pearsonr(sub["Discharge_mean"], sub[col])
    corr_results.append({
        "Supply_Type": col,
        "Pearson_r":   round(r, 4),
        "P_value":     round(p, 4),
        "N_HUC4":      len(sub),
        "Significant": "Yes" if p < 0.05 else "No"
    })
    print(f"  {col:15s} | r = {r:.4f} | p = {p:.4f} | n = {len(sub)} | {'✓ Significant' if p < 0.05 else '✗ Not significant'}")

corr_df = pd.DataFrame(corr_results)

# ============================================================
# STEP 7: Save CSV + Shapefile
# ============================================================
print("\n" + "=" * 50)
print("STEP 7: Saving results...")
print("=" * 50)

# Save correlation summary
corr_csv = os.path.join(OUTPUT_DIR, "correlation_discharge_supply.csv")
corr_df.to_csv(corr_csv, index=False)
print(f"✓ Correlation summary: {corr_csv}")

# Save full HUC4 dataset
full_csv = os.path.join(OUTPUT_DIR, "HUC4_discharge_supply.csv")
merged.to_csv(full_csv, index=False)
print(f"✓ Full dataset: {full_csv}")

# Save shapefile (join with HUC4 geometry)
shp_merged = huc4_gdf.merge(
    merged,
    left_on  = "huc4",
    right_on = "HUC4",
    how="inner"
)
shp_out = os.path.join(OUTPUT_DIR, "HUC4_discharge_supply.shp")
shp_merged.to_file(shp_out)
print(f"✓ Shapefile: {shp_out}")

print("\n" + "=" * 50)
print("All done! Results in output/correlation/")
print("=" * 50)
