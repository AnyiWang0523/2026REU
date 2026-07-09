# ============================================================
# Monthly Supply Data Export
# Category: Public Supply - Surface Water
# Output: HUC4 | Year | Month | SW_Supply (Mgal/month)
# ============================================================

import pandas as pd
import os

# ============================================================
# Path Configuration
# ============================================================

BASE_DIR   = r"C:\Users\anyiw\OneDrive - University of Illinois - Urbana\GIS\2026REU\WSDI"
CSV_PATH   = os.path.join(BASE_DIR, "data", "CSV", "PS_HUC12_SW_2000_2020.csv")
OUTPUT_DIR = os.path.join(BASE_DIR, "output", "PS_SW", "monthly")
os.makedirs(OUTPUT_DIR, exist_ok=True)

# ============================================================
# STEP 1: Load Data
# ============================================================
print("=" * 50)
print("STEP 1: Loading PS SW data from local CSV...")
print("=" * 50)

sw_raw = pd.read_csv(CSV_PATH, dtype=str)
print(f"✓ Load complete: {sw_raw.shape[0]} rows x {sw_raw.shape[1]} columns")

# ============================================================
# STEP 2: Parse CSV Structure
# ============================================================
print("\n" + "=" * 50)
print("STEP 2: Parsing CSV structure...")
print("=" * 50)

all_cols   = sw_raw.columns.tolist()
meta_cols  = ["Year", "Month"]
huc12_cols = all_cols[2:]

def fix_huc12(val):
    """Convert scientific notation column names to 12-digit HUC12 strings"""
    try:
        return str(int(float(str(val).strip()))).zfill(12)
    except:
        return str(val).strip().zfill(12)

huc12_fixed    = [fix_huc12(c) for c in huc12_cols]
sw_raw.columns = meta_cols + huc12_fixed
sw_raw["Year"]  = sw_raw["Year"].astype(int)
sw_raw["Month"] = sw_raw["Month"].astype(int)

print(f"Number of HUC12 units: {len(huc12_fixed)}")

# ============================================================
# STEP 3: Wide to Long Format -> Aggregate to HUC4 Monthly
# ============================================================
print("\n" + "=" * 50)
print("STEP 3: Aggregating to HUC4 monthly totals...")
print("=" * 50)

sw_long = sw_raw.melt(
    id_vars    = ["Year", "Month"],
    value_vars = huc12_fixed,
    var_name   = "HUC12",
    value_name = "SW_Mgal_day"
)

sw_long["SW_Mgal_day"] = pd.to_numeric(
    sw_long["SW_Mgal_day"], errors="coerce"
).fillna(0)

sw_long = sw_long[sw_long["SW_Mgal_day"] > 0].copy()

# Convert Mgal/day -> Mgal/month
days_map = {1:31, 2:28, 3:31, 4:30, 5:31,  6:30,
            7:31, 8:31, 9:30, 10:31, 11:30, 12:31}
sw_long["SW_Mgal_month"] = sw_long["SW_Mgal_day"] * sw_long["Month"].map(days_map)

# Extract HUC4 from first 4 digits of HUC12
sw_long["HUC4"] = sw_long["HUC12"].str[:4]

# Sum to HUC4 + Year + Month
sw_monthly = (
    sw_long
    .groupby(["HUC4", "Year", "Month"])["SW_Mgal_month"]
    .sum()
    .reset_index()
    .rename(columns={"SW_Mgal_month": "SW_Supply_Mgal"})
)

print(f"✓ Monthly aggregation complete: {len(sw_monthly)} records")
print(f"  Unique HUC4 units: {sw_monthly['HUC4'].nunique()}")
print(f"  Year range: {sw_monthly['Year'].min()} - {sw_monthly['Year'].max()}")
print(sw_monthly.head(10).to_string(index=False))

# ============================================================
# STEP 4: Save Results
# ============================================================
print("\n" + "=" * 50)
print("STEP 4: Saving results...")
print("=" * 50)

# Save one CSV per year
for yr in range(2000, 2021):
    sub = sw_monthly[sw_monthly["Year"] == yr].copy()
    if sub.empty:
        continue
    out_csv = os.path.join(OUTPUT_DIR, f"PS_SW_monthly_{yr}.csv")
    sub.to_csv(out_csv, index=False)
    print(f"  ✓ PS_SW_monthly_{yr}.csv saved ({len(sub)} records)")

# Also save complete dataset as single CSV
out_full = os.path.join(OUTPUT_DIR, "PS_SW_monthly_2000_2020.csv")
sw_monthly.to_csv(out_full, index=False)
print(f"\n✓ Full dataset saved: PS_SW_monthly_2000_2020.csv")

print("\n" + "=" * 50)
print("All done! Results are in output/PS_SW/monthly/ folder")
print("=" * 50)