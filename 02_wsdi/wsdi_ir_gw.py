# ============================================================
# WSDI (Water Supply-Demand Index) Calculation
# Category: Irrigation - Groundwater
# Formula: (GW_Supply - Demand) / Demand
# Supply: IR_HUC12_GW_WD_m_2000_2020.csv (local CSV)
# Demand: Supply value from the year when PDSI is closest to 0
# Spatial unit: HUC4
# ============================================================

import geopandas as gpd
import pandas as pd
import numpy as np
import requests
import time
import os

# ============================================================
# Path Configuration
# ============================================================

BASE_DIR   = r"C:\Users\anyiw\OneDrive - University of Illinois - Urbana\GIS\2026REU\WSDI"
HUC4_SHP   = os.path.join(BASE_DIR, "data", "HUC4",   "WBDHU4.shp")
STATE_SHP  = os.path.join(BASE_DIR, "data", "States", "cb_2020_us_state_500k.shp")
CSV_PATH   = os.path.join(BASE_DIR, "data", "CSV",    "IR_HUC12_GW_WD_monthly_2000_2020.csv")
OUTPUT_DIR = os.path.join(BASE_DIR, "output", "IR_GW")
os.makedirs(OUTPUT_DIR, exist_ok=True)

# ============================================================
# STEP 1: Load GW Data from Local CSV
# ============================================================
print("=" * 50)
print("STEP 1: Loading IR GW data from local CSV...")
print("=" * 50)

gw_raw = pd.read_csv(CSV_PATH, dtype=str)
print(f"✓ Load complete: {gw_raw.shape[0]} rows x {gw_raw.shape[1]} columns")
print(f"First 5 columns: {gw_raw.columns.tolist()[:5]}")
print(f"Sample rows:")
print(gw_raw.head(3).to_string())

# ============================================================
# STEP 2: Parse CSV Structure
# ============================================================
print("\n" + "=" * 50)
print("STEP 2: Parsing CSV structure...")
print("=" * 50)

all_cols   = gw_raw.columns.tolist()
meta_cols  = ["Year", "Month"]
huc12_cols = all_cols[2:]

def fix_huc12(val):
    """Convert scientific notation column names to 12-digit HUC12 strings"""
    try:
        return str(int(float(str(val).strip()))).zfill(12)
    except:
        return str(val).strip().zfill(12)

huc12_fixed    = [fix_huc12(c) for c in huc12_cols]
gw_raw.columns = meta_cols + huc12_fixed
gw_raw["Year"]  = gw_raw["Year"].astype(int)
gw_raw["Month"] = gw_raw["Month"].astype(int)

print(f"Number of HUC12 units: {len(huc12_fixed)}")
print(f"HUC12 sample: {huc12_fixed[:3]}")

# ============================================================
# STEP 3: Wide to Long Format -> Aggregate to HUC4 Annual
# ============================================================
print("\n" + "=" * 50)
print("STEP 3: Aggregating to HUC4 annual totals...")
print("=" * 50)

gw_long = gw_raw.melt(
    id_vars    = ["Year", "Month"],
    value_vars = huc12_fixed,
    var_name   = "HUC12",
    value_name = "GW_Mgal_day"
)

gw_long["GW_Mgal_day"] = pd.to_numeric(
    gw_long["GW_Mgal_day"], errors="coerce"
).fillna(0)

gw_long = gw_long[gw_long["GW_Mgal_day"] > 0].copy()

days_map = {1:31, 2:28, 3:31, 4:30, 5:31,  6:30,
            7:31, 8:31, 9:30, 10:31, 11:30, 12:31}
gw_long["GW_Mgal"] = gw_long["GW_Mgal_day"] * gw_long["Month"].map(days_map)
gw_long["HUC4"]    = gw_long["HUC12"].str[:4]

gw_annual = (
    gw_long
    .groupby(["HUC4", "Year"])["GW_Mgal"]
    .sum()
    .reset_index()
    .rename(columns={"GW_Mgal": "GW_Supply"})
)

print(f"✓ Annual aggregation complete: {len(gw_annual)} records")
print(f"  Unique HUC4 units: {gw_annual['HUC4'].nunique()}")
print(gw_annual.head(5).to_string(index=False))

# ============================================================
# STEP 4: Spatial Join HUC4 -> State
# ============================================================
print("\n" + "=" * 50)
print("STEP 4: Spatial join HUC4 -> State...")
print("=" * 50)

huc4_gdf  = gpd.read_file(HUC4_SHP)
state_gdf = gpd.read_file(STATE_SHP)
state_gdf = state_gdf.to_crs(huc4_gdf.crs)

huc4_col  = "huc4"
state_col = "STUSPS"

joined = gpd.sjoin(
    huc4_gdf[[huc4_col, "geometry"]],
    state_gdf[[state_col, "geometry"]],
    how="left",
    predicate="intersects"
)

huc4_state = (
    joined[[huc4_col, state_col]]
    .dropna()
    .drop_duplicates(subset=huc4_col, keep="first")
    .rename(columns={huc4_col: "HUC4", state_col: "STATE"})
)
huc4_state["HUC4"] = huc4_state["HUC4"].astype(str).str.zfill(4)

print(f"✓ HUC4-State mapping complete: {len(huc4_state)} records")
print(huc4_state.head(5).to_string(index=False))

# ============================================================
# STEP 5: Fetch Annual PDSI by State (NOAA API)
# ============================================================
print("\n" + "=" * 50)
print("STEP 5: Fetching PDSI data from NOAA API...")
print("=" * 50)

NOAA_CODES = {
    "AL":"01", "AZ":"02", "AR":"03", "CA":"04", "CO":"05", "CT":"06",
    "DE":"07", "FL":"08", "GA":"09", "ID":"10", "IL":"11", "IN":"12",
    "IA":"13", "KS":"14", "KY":"15", "LA":"16", "ME":"17", "MD":"18",
    "MA":"19", "MI":"20", "MN":"21", "MS":"22", "MO":"23", "MT":"24",
    "NE":"25", "NV":"26", "NH":"27", "NJ":"28", "NM":"29", "NY":"30",
    "NC":"31", "ND":"32", "OH":"33", "OK":"34", "OR":"35", "PA":"36",
    "RI":"37", "SC":"38", "SD":"39", "TN":"40", "TX":"41", "UT":"42",
    "VT":"43", "VA":"44", "WA":"45", "WV":"46", "WI":"47", "WY":"48"
}

def fetch_pdsi(state_abbr, beg=2000, end=2020):
    """Fetch annual PDSI from NOAA Climate at a Glance API"""
    code = NOAA_CODES.get(state_abbr.upper())
    if not code:
        return pd.DataFrame()
    url = (
        f"https://www.ncei.noaa.gov/access/monitoring/climate-at-a-glance/"
        f"statewide/time-series/{code}/pdsi/all/0/{beg}-{end}/data.json"
    )
    try:
        r = requests.get(url, timeout=30)
        r.raise_for_status()
        rows = [
            {"STATE": state_abbr.upper(),
             "Year":  int(k[:4]),
             "PDSI":  float(v["value"])}
            for k, v in r.json()["data"].items()
        ]
        return pd.DataFrame(rows)
    except Exception as e:
        print(f"  ✗ {state_abbr}: {e}")
        return pd.DataFrame()

states_needed = huc4_state["STATE"].dropna().unique().tolist()
print(f"Fetching PDSI for {len(states_needed)} states...")

pdsi_list = []
for st in states_needed:
    df = fetch_pdsi(st)
    if not df.empty:
        pdsi_list.append(df)
        best = df.loc[df["PDSI"].abs().idxmin()]
        print(f"  ✓ {st} | Reference year: {int(best['Year'])} | PDSI = {best['PDSI']:.3f}")
    time.sleep(0.3)

pdsi_all = pd.concat(pdsi_list, ignore_index=True)
print(f"✓ PDSI fetch complete: {len(pdsi_all)} records")

# ============================================================
# STEP 6: Find Demand Reference Year + Calculate WSDI
# ============================================================
print("\n" + "=" * 50)
print("STEP 6: Calculating WSDI...")
print("=" * 50)

demand_ref = (
    pdsi_all
    .assign(abs_pdsi=lambda x: x["PDSI"].abs())
    .sort_values("abs_pdsi")
    .drop_duplicates("STATE")
    [["STATE", "Year", "PDSI"]]
    .rename(columns={"Year": "Demand_Year", "PDSI": "Demand_PDSI"})
)

print("\nDemand reference year per state (year with PDSI closest to 0):")
print(demand_ref.sort_values("STATE").to_string(index=False))

huc4_ref = huc4_state.merge(demand_ref, on="STATE", how="left")

huc4_demand = huc4_ref.merge(
    gw_annual,
    left_on  = ["HUC4", "Demand_Year"],
    right_on = ["HUC4", "Year"],
    how="left"
).rename(columns={"GW_Supply": "Demand"})[
    ["HUC4", "STATE", "Demand_Year", "Demand_PDSI", "Demand"]
]

wsdi = (
    gw_annual
    .merge(huc4_demand, on="HUC4", how="left")
    .query("Demand > 0 and Demand == Demand")
    .copy()
)

wsdi["WSDI"] = (wsdi["GW_Supply"] - wsdi["Demand"]) / wsdi["Demand"]

print(f"\n✓ WSDI calculation complete: {len(wsdi)} records")
print("\nWSDI summary statistics:")
print(wsdi["WSDI"].describe().round(4))
print("\nSample results:")
print(wsdi[["HUC4","STATE","Year","GW_Supply","Demand","Demand_Year","WSDI"]].head(10).to_string(index=False))

# ============================================================
# STEP 7: Save Results
# ============================================================
print("\n" + "=" * 50)
print("STEP 7: Saving results...")
print("=" * 50)

out_csv = os.path.join(OUTPUT_DIR, "WSDI_IR_GW_HUC4_2000_2020.csv")
wsdi.to_csv(out_csv, index=False)
print(f"✓ Full CSV saved: {out_csv}")

for yr in range(2000, 2021):
    sub = wsdi[wsdi["Year"] == yr]
    if sub.empty:
        continue

    merged = huc4_gdf.merge(
        sub[["HUC4", "WSDI", "GW_Supply", "Demand", "Demand_Year", "STATE"]],
        left_on  = huc4_col,
        right_on = "HUC4",
        how="left"
    )

    out_shp = os.path.join(OUTPUT_DIR, f"HUC4_WSDI_IR_GW_{yr}.shp")
    merged.to_file(out_shp)
    print(f"  ✓ HUC4_WSDI_IR_GW_{yr}.shp saved")

print("\n" + "=" * 50)
print("All done! Results are in output/IR_GW/ folder")
print("=" * 50)