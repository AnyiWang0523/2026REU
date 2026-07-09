# ============================================================
# Monthly Supply Data Export
# Category: Public Supply - Groundwater
# Output: HUC4 | Year | Month | GW_Supply (Mgal/month)
# ============================================================

import pandas as pd
import os

# ============================================================
# Path Configuration
# ============================================================

BASE_DIR   = r"C:\Users\anyiw\OneDrive - University of Illinois - Urbana\GIS\2026REU\WSDI"
CSV_PATH   = os.path.join(BASE_DIR, "data", "CSV", "PS_HUC12_GW_2000_2020.csv")
OUTPUT_DIR = os.path.join(BASE_DIR, "output", "PS_GW", "monthly")
os.makedirs(OUTPUT_DIR, exist_ok=True)

# ============================================================
# STEP 1: Load Data
# ============================================================
print("=" * 50)
print("STEP 1: Loading PS GW data from local CSV...")
print("=" * 50)

gw_raw = pd.read_csv(CSV_PATH, dtype=str)
print(f"✓ Load complete: {gw_raw.shape[0]} rows x {gw_raw.shape[1]} columns")

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

# ============================================================
# STEP 3: Wide to Long Format -> Aggregate to HUC4 Monthly
# ============================================================
print("\n" + "=" * 50)
print("STEP 3: Aggregating to HUC4 monthly totals...")
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

# Convert Mgal/day -> Mgal/month
days_map = {1:31, 2:28, 3:31, 4:30, 5:31,  6:30,
            7:31, 8:31, 9:30, 10:31, 11:30, 12:31}
gw_long["GW_Mgal_month"] = gw_long["GW_Mgal_day"] * gw_long["Month"].map(days_map)

# Extract HUC4 from first 4 digits of HUC12
gw_long["HUC4"] = gw_long["HUC12"].str[:4]

# Sum to HUC4 + Year + Month
gw_monthly = (
    gw_long
    .groupby(["HUC4", "Year", "Month"])["GW_Mgal_month"]
    .sum()
    .reset_index()
    .rename(columns={"GW_Mgal_month": "GW_Supply_Mgal"})
)

print(f"✓ Monthly aggregation complete: {len(gw_monthly)} records")
print(f"  Unique HUC4 units: {gw_monthly['HUC4'].nunique()}")
print(f"  Year range: {gw_monthly['Year'].min()} - {gw_monthly['Year'].max()}")
print(gw_monthly.head(10).to_string(index=False))

# ============================================================
# STEP 4: Save Results
# ============================================================
print("\n" + "=" * 50)
print("STEP 4: Saving results...")
print("=" * 50)

# Save one CSV per year
for yr in range(2000, 2021):
    sub = gw_monthly[gw_monthly["Year"] == yr].copy()
    if sub.empty:
        continue
    out_csv = os.path.join(OUTPUT_DIR, f"PS_GW_monthly_{yr}.csv")
    sub.to_csv(out_csv, index=False)
    print(f"  ✓ PS_GW_monthly_{yr}.csv saved ({len(sub)} records)")

# Also save complete dataset as single CSV
out_full = os.path.join(OUTPUT_DIR, "PS_GW_monthly_2000_2020.csv")
gw_monthly.to_csv(out_full, index=False)
print(f"\n✓ Full dataset saved: PS_GW_monthly_2000_2020.csv")

print("\n" + "=" * 50)
print("All done! Results are in output/PS_GW/monthly/ folder")
print("=" * 50)
