import os
import pandas as pd

BASE_DIR   = r"C:\Users\anyiw\OneDrive - University of Illinois - Urbana\GIS\2026REU\WSDI"
COV_CSV    = os.path.join(BASE_DIR, "output", "reservoir_coverage_table.csv")
FLAG_CSV   = os.path.join(BASE_DIR, "output", "reservoir_storage_flags.csv")
OUTPUT     = os.path.join(BASE_DIR, "output", "huc4_valid_coverage50.csv")

PCT_COL    = "Col5_Coverage_Pct_NIDMatched_vs_NIDTotal"
PCT_THRESH = 50.0

# ============================================================
# Load
# ============================================================
cov   = pd.read_csv(COV_CSV, dtype={"HUC4": str})
flags = pd.read_csv(FLAG_CSV, dtype={"huc4": str})

print(f"Coverage table rows : {len(cov)}")
print(f"Flag rows           : {len(flags)}")

# ============================================================
# Step 1 — remove any HUC4 that appears in the flags file
# ============================================================
flagged_huc4s = set(flags["huc4"].dropna().unique())
print(f"\nDistinct flagged HUC4s : {len(flagged_huc4s)}")

cov_clean = cov[~cov["HUC4"].isin(flagged_huc4s)].copy()
print(f"After removing flagged : {len(cov_clean)} rows")

# ============================================================
# Step 2 — keep only HUC4s with 50% <= Col5b <= 100%
#   - >= 50%  : GRanD covers at least half of matched NID capacity
#   - <= 100% : safety net (coverage_over100 should already be in flags,
#               but guard against re-runs with stale flag files)
# ============================================================
cov_clean[PCT_COL] = pd.to_numeric(cov_clean[PCT_COL], errors="coerce")
cov_valid = cov_clean[
    (cov_clean[PCT_COL] >= PCT_THRESH) &
    (cov_clean[PCT_COL] <= 100.0)
].copy()
print(f"After 50% <= Col5b <= 100% : {len(cov_valid)} rows")

# ============================================================
# Save
# ============================================================
cov_valid.to_csv(OUTPUT, index=False, encoding="utf-8-sig")
print(f"\nSaved -> {OUTPUT}")

print("\n--- Preview ---")
print(cov_valid[["HUC4", "HUC4_Name", PCT_COL]].to_string(index=False))
