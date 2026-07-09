# ============================================================
# Fetch Annual Riverflow for HUC4s with Reservoirs
# Only HUC4s with large reservoirs (missing_ratio < 0.05)
# Output: annual mean streamflow per HUC4 (2000-2020)
# ============================================================

import dataretrieval.nwis as nwis
import geopandas as gpd
import pandas as pd
import numpy as np
import os
import time
from scipy.spatial import cKDTree

# ============================================================
# Path Configuration
# ============================================================

BASE_DIR   = r"C:\Users\anyiw\OneDrive - University of Illinois - Urbana\GIS\2026REU\WSDI"
HUC4_SHP   = os.path.join(BASE_DIR, "data", "HUC4", "WBDHU4.shp")
RESV_CSV   = os.path.join(BASE_DIR, "data", "CSV", "GDROM_Resv.csv")
OUTPUT_DIR = os.path.join(BASE_DIR, "output", "correlation")
os.makedirs(OUTPUT_DIR, exist_ok=True)

# ============================================================
# STEP 1: Load HUC4 boundaries
# ============================================================
print("=" * 50)
print("STEP 1: Loading HUC4 boundaries...")
print("=" * 50)

huc4_gdf = gpd.read_file(HUC4_SHP)
huc4_gdf["huc4"] = huc4_gdf["huc4"].astype(str).str.zfill(4)
huc4_gdf = huc4_gdf.to_crs("EPSG:4326")
huc4_gdf["centroid_lon"] = huc4_gdf.geometry.centroid.x
huc4_gdf["centroid_lat"] = huc4_gdf.geometry.centroid.y
huc4_list = huc4_gdf["huc4"].tolist()
print(f"✓ HUC4 units loaded: {len(huc4_list)}")

# ============================================================
# STEP 2: Find HUC4s with reservoirs
# ============================================================
print("\n" + "=" * 50)
print("STEP 2: Finding HUC4s with large reservoirs...")
print("=" * 50)

resv = pd.read_csv(RESV_CSV)
resv["missing_ratio"] = pd.to_numeric(resv["missing_ratio"], errors="coerce")
resv = resv[resv["missing_ratio"] < 0.05].copy()
resv = resv.dropna(subset=["lattitude", "longtitude"])

resv_gdf = gpd.GeoDataFrame(
    resv,
    geometry=gpd.points_from_xy(resv["longtitude"], resv["lattitude"]),
    crs="EPSG:4326"
)

resv_joined = gpd.sjoin(
    resv_gdf,
    huc4_gdf[["huc4", "geometry"]],
    how="left",
    predicate="within"
)

huc4_with_resv = sorted(resv_joined["huc4"].dropna().unique().tolist())
print(f"✓ HUC4s with reservoirs: {len(huc4_with_resv)}")
print(f"  HUC4 list: {huc4_with_resv[:10]}...")

# ============================================================
# STEP 3: Fetch streamflow stations for these HUC4s only
# ============================================================
print("\n" + "=" * 50)
print("STEP 3: Fetching streamflow stations...")
print("=" * 50)

huc2_list = sorted(set([h[:2] for h in huc4_with_resv]))
print(f"  Querying {len(huc2_list)} HUC2 regions...")

all_station_data = []

for huc2 in huc2_list:
    try:
        sites = nwis.get_info(
            huc=huc2,
            parameterCd="00060",
            siteType="ST"
        )[0]

        if sites.empty:
            continue

        sites["HUC4"] = sites["huc_cd"].astype(str).str.zfill(8).str[:4]

        # Only keep stations in HUC4s with reservoirs
        sites = sites[sites["HUC4"].isin(huc4_with_resv)]

        if sites.empty:
            continue

        all_station_data.append(sites[["site_no", "station_nm",
                                        "dec_lat_va", "dec_long_va", "HUC4"]])
        print(f"  ✓ HUC2 {huc2}: {len(sites)} stations across {sites['HUC4'].nunique()} HUC4s")

    except Exception as e:
        print(f"  ✗ HUC2 {huc2}: {e}")

    time.sleep(0.3)

if not all_station_data:
    print("✗ No stations found!")
    exit()

stations_df = pd.concat(all_station_data, ignore_index=True)
print(f"\n✓ Total stations before filter: {len(stations_df)}")

# Max 3 stations per HUC4
stations_df = (
    stations_df
    .groupby("HUC4")
    .head(3)
    .reset_index(drop=True)
)
print(f"✓ Stations after filter (max 3/HUC4): {len(stations_df)}")
print(f"  HUC4s with stations: {stations_df['HUC4'].nunique()}")

stations_df.to_csv(os.path.join(OUTPUT_DIR, "resv_streamflow_stations.csv"), index=False)

# ============================================================
# STEP 4: Fetch annual streamflow 2000-2020
# ============================================================
print("\n" + "=" * 50)
print("STEP 4: Fetching annual streamflow data...")
print("=" * 50)

annual_records = []
total_stations = len(stations_df)

for count, (i, row) in enumerate(stations_df.iterrows()):
    site = str(row["site_no"])
    huc4 = row["HUC4"]
    try:
        df, _ = nwis.get_dv(
            sites=site,
            parameterCd="00060",
            start="2000-01-01",
            end="2020-12-31"
        )

        if df.empty:
            continue

        flow_col = [c for c in df.columns if "00060" in c and "Mean" in c]
        if not flow_col:
            flow_col = [c for c in df.columns if "00060" in c]
        if not flow_col:
            continue

        df = df[[flow_col[0]]].rename(columns={flow_col[0]: "flow_cfs"})
        df["flow_cfs"] = pd.to_numeric(df["flow_cfs"], errors="coerce")
        df = df.dropna()
        df.index = pd.to_datetime(df.index)
        df["Year"] = df.index.year

        annual = df.groupby("Year")["flow_cfs"].mean().reset_index()
        annual = annual[(annual["Year"] >= 2000) & (annual["Year"] <= 2020)]
        annual["site_no"] = site
        annual["HUC4"]    = huc4
        annual_records.append(annual)

    except Exception as e:
        pass

    time.sleep(0.1)

    if (count + 1) % 10 == 0:
        print(f"  Progress: {count+1}/{total_stations} stations fetched")

if not annual_records:
    print("✗ No flow data retrieved!")
    exit()

flow_df = pd.concat(annual_records, ignore_index=True)
print(f"✓ Annual flow records: {len(flow_df)}")

# ============================================================
# STEP 5: Aggregate to HUC4 annual mean
# ============================================================
print("\n" + "=" * 50)
print("STEP 5: Aggregating to HUC4 annual mean...")
print("=" * 50)

huc4_flow = (
    flow_df
    .groupby(["HUC4", "Year"])["flow_cfs"]
    .mean()
    .reset_index()
    .rename(columns={"flow_cfs": "Flow_cfs"})
)

huc4_flow["Flow_Mgal_day"] = huc4_flow["Flow_cfs"] * 0.646317

print(f"✓ HUC4 annual flow: {len(huc4_flow)} records")
print(f"  HUC4s with data: {huc4_flow['HUC4'].nunique()}")

# ============================================================
# STEP 6: Correlation with PS GW and SW Supply
# ============================================================
print("\n" + "=" * 50)
print("STEP 6: Correlation Analysis...")
print("=" * 50)

from scipy import stats

PS_GW_CSV = os.path.join(BASE_DIR, "output", "PS_GW", "WSDI_PS_GW_HUC4_2000_2020.csv")
PS_SW_CSV = os.path.join(BASE_DIR, "output", "PS_SW", "WSDI_PS_SW_HUC4_2000_2020.csv")

ps_gw = pd.read_csv(PS_GW_CSV, dtype={"HUC4": str})
ps_sw = pd.read_csv(PS_SW_CSV, dtype={"HUC4": str})
ps_gw["HUC4"] = ps_gw["HUC4"].str.zfill(4)
ps_sw["HUC4"] = ps_sw["HUC4"].str.zfill(4)

# Keep only HUC4s with reservoirs
ps_gw = ps_gw[ps_gw["HUC4"].isin(huc4_with_resv)]
ps_sw = ps_sw[ps_sw["HUC4"].isin(huc4_with_resv)]

# Merge flow with supply by HUC4 + Year
merged_gw = huc4_flow.merge(
    ps_gw[["HUC4", "Year", "GW_Supply"]],
    on=["HUC4", "Year"], how="inner"
)
merged_sw = huc4_flow.merge(
    ps_sw[["HUC4", "Year", "SW_Supply"]],
    on=["HUC4", "Year"], how="inner"
)

print(f"  GW merged records: {len(merged_gw)}")
print(f"  SW merged records: {len(merged_sw)}")

# Overall correlation
r_gw, p_gw = stats.pearsonr(merged_gw["Flow_cfs"], merged_gw["GW_Supply"])
r_sw, p_sw = stats.pearsonr(merged_sw["Flow_cfs"], merged_sw["SW_Supply"])

print(f"\n  Discharge vs PS GW | r = {r_gw:.4f} | p = {p_gw:.4f} | {'✓ Significant' if p_gw < 0.05 else '✗ Not significant'}")
print(f"  Discharge vs PS SW | r = {r_sw:.4f} | p = {p_sw:.4f} | {'✓ Significant' if p_sw < 0.05 else '✗ Not significant'}")

# Per HUC4 correlation
huc4_corr = []
for huc4 in huc4_with_resv:
    sub_gw = merged_gw[merged_gw["HUC4"] == huc4]
    sub_sw = merged_sw[merged_sw["HUC4"] == huc4]

    row = {"HUC4": huc4}

    if len(sub_gw) >= 5:
        r, p = stats.pearsonr(sub_gw["Flow_cfs"], sub_gw["GW_Supply"])
        row["r_GW"] = round(r, 4)
        row["p_GW"] = round(p, 4)
    else:
        row["r_GW"] = np.nan
        row["p_GW"] = np.nan

    if len(sub_sw) >= 5:
        r, p = stats.pearsonr(sub_sw["Flow_cfs"], sub_sw["SW_Supply"])
        row["r_SW"] = round(r, 4)
        row["p_SW"] = round(p, 4)
    else:
        row["r_SW"] = np.nan
        row["p_SW"] = np.nan

    huc4_corr.append(row)

corr_df = pd.DataFrame(huc4_corr)
corr_df["sig_GW"] = corr_df["p_GW"].apply(lambda x: "Yes" if x < 0.05 else "No")
corr_df["sig_SW"] = corr_df["p_SW"].apply(lambda x: "Yes" if x < 0.05 else "No")

print(f"\n  HUC4s with significant GW correlation: {(corr_df['sig_GW']=='Yes').sum()}")
print(f"  HUC4s with significant SW correlation: {(corr_df['sig_SW']=='Yes').sum()}")

# ============================================================
# STEP 7: Save Results
# ============================================================
print("\n" + "=" * 50)
print("STEP 7: Saving results...")
print("=" * 50)

# Save flow data
flow_csv = os.path.join(OUTPUT_DIR, "riverflow_resv_HUC4_annual_2000_2020.csv")
huc4_flow.to_csv(flow_csv, index=False)
print(f"✓ Flow data: {flow_csv}")

# Save correlation results
corr_csv = os.path.join(OUTPUT_DIR, "correlation_discharge_supply_annual.csv")
corr_df.to_csv(corr_csv, index=False)
print(f"✓ Correlation results: {corr_csv}")

# Save shapefile
shp_merged = huc4_gdf.merge(corr_df, left_on="huc4", right_on="HUC4", how="inner")
shp_out = os.path.join(OUTPUT_DIR, "HUC4_discharge_supply_corr.shp")
shp_merged.to_file(shp_out)
print(f"✓ Shapefile: {shp_out}")

print("\n" + "=" * 50)
print("All done! Results in output/correlation/")
print("=" * 50)