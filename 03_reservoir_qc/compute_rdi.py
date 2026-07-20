"""
compute_rdi.py  -  RDI (Reservoir Drought Index), all valid HUC4s
=================================================================
Implements the R-code formula from Vora & Cai (2026):

    Drought_Threshold[m] = quantile(monthly_avg_storage, 0.25)   per calendar month
    RDI_t = (monthly_avg_storage_t - Drought_Threshold[month(t)]) / mean(all monthly_avg)

HUC4 aggregation by aggregate-first approach:
    S_HUC4_t    = sum_i(S_i_t)
    RDI_HUC4_t  = (S_HUC4_t - Threshold[month(t)]) / mean(S_HUC4)

Data source:
    data/hydroshare_data/{GRanD_ID}.csv
    Columns: Time(YYYY-MM-DD), Storage(fraction 0-1), NetInflow(fraction), Release(fraction)
    Reservoir IDs and capacities from data/CSV/GDROM_Resv.csv

Analysis period: dynamic — uses only months where ALL reservoirs in a HUC4 have valid data.

Outputs:
    output/rdi/all_HUC4_RDI.csv            -- monthly RDI for all valid HUC4s
    output/rdi/all_reservoir_monthly.csv    -- per-reservoir monthly storage + RDI
"""

import os
import sys
import warnings
import numpy as np
import pandas as pd
import geopandas as gpd

warnings.filterwarnings("ignore")
sys.stdout.reconfigure(encoding="utf-8", errors="replace")

# ── Paths ─────────────────────────────────────────────────────────────────────
BASE_DIR       = r"C:\Users\anyiw\OneDrive - University of Illinois - Urbana\GIS\2026REU\WSDI"
GDROM_CSV      = os.path.join(BASE_DIR, "data", "CSV", "GDROM_Resv.csv")
HUC4_SHP       = os.path.join(BASE_DIR, "data", "HUC4", "WBDHU4.shp")
HYDROSHARE_DIR = os.path.join(BASE_DIR, "data", "hydroshare_data")
VALID_CSV      = os.path.join(BASE_DIR, "output", "huc4_valid_coverage50.csv")
OUTPUT_DIR     = os.path.join(BASE_DIR, "output", "rdi")
os.makedirs(OUTPUT_DIR, exist_ok=True)

AF_TO_CFS = 43560 / 86400

# ── Load metadata once ────────────────────────────────────────────────────────
print("=" * 60)
print("Loading metadata ...")
print("=" * 60)

gdrom = pd.read_csv(GDROM_CSV).rename(columns={
    "lattitude": "latitude", "longtitude": "longitude"
})
hs_ids = {int(f.replace(".csv", "")) for f in os.listdir(HYDROSHARE_DIR)
          if f.endswith(".csv")}
gdrom["has_hs"] = gdrom["GRanD_ID"].isin(hs_ids)

gdrom_gdf = gpd.GeoDataFrame(
    gdrom, crs="EPSG:4326",
    geometry=gpd.points_from_xy(gdrom["longitude"], gdrom["latitude"])
)

huc4_gdf = gpd.read_file(HUC4_SHP).to_crs("EPSG:4326")
huc4_gdf["huc4"] = huc4_gdf["huc4"].astype(str).str.zfill(4)

valid_huc4 = pd.read_csv(VALID_CSV, dtype={"HUC4": str})
valid_huc4["HUC4"] = valid_huc4["HUC4"].str.zfill(4)

# Spatial join: assign each GDROM reservoir to its HUC4
all_joined = gpd.sjoin(
    gdrom_gdf, huc4_gdf[["huc4", "geometry"]],
    how="left", predicate="within"
).reset_index(drop=True)

print(f"  Valid HUC4s          : {len(valid_huc4)}")
print(f"  HydroShare files     : {len(hs_ids)}")
print(f"  GDROM with HS data   : {gdrom['has_hs'].sum()}")

# ── Core functions ────────────────────────────────────────────────────────────

MAX_MISSING_FRAC = 0.2   # months with >20% missing days are excluded from RDI

def load_hydroshare(grand_id: int, cap_af: float) -> pd.Series:
    """
    Read {grand_id}.csv, convert Storage fraction -> acre-feet,
    interpolate missing values, then mask months that have too many
    missing days as NaN so they are excluded from RDI calculation.

    A month is kept only if: missing_days / total_days <= MAX_MISSING_FRAC.
    Interpolation still fills individual gaps for plotting continuity,
    but the monthly average is set to NaN for incomplete months.
    """
    path = os.path.join(HYDROSHARE_DIR, f"{grand_id}.csv")
    raw  = pd.read_csv(path, usecols=["Time", "Storage"])
    raw["Time"] = pd.to_datetime(raw["Time"])
    raw = raw.drop_duplicates(subset="Time").set_index("Time").sort_index()

    dates = pd.date_range(raw.index[0], raw.index[-1], freq="D")
    strg_raw = raw["Storage"].reindex(dates) * cap_af   # fraction -> af; NaN where missing

    # Per-month completeness: fraction of days with real data
    month_idx       = strg_raw.resample("MS").count()       # real days count
    month_total     = strg_raw.resample("MS").size()        # total days in month
    month_miss_frac = 1 - month_idx / month_total           # missing fraction

    # Months that are too incomplete -> will be NaN after resampling
    bad_months = set(month_miss_frac[month_miss_frac > MAX_MISSING_FRAC].index)

    # Interpolate for daily continuity, then re-mask bad months
    strg = strg_raw.clip(lower=0).interpolate(method="linear", limit_direction="both")

    # Zero out storage in bad months so monthly mean becomes NaN-like;
    # easier: set bad-month days back to NaN after interpolation
    for bm in bad_months:
        mask = (strg.index.year == bm.year) & (strg.index.month == bm.month)
        strg[mask] = np.nan

    return strg


def compute_rdi(monthly_s: pd.Series) -> tuple[pd.Series, np.ndarray]:
    """
    R-code formula:
        Threshold[m] = quantile(S where Month==m, 0.25)
        RDI_t = (S_t - Threshold[month(t)]) / mean(S_all)
    Returns (RDI Series on monthly_s.index, thresholds[12])
    """
    overall_mean = monthly_s.mean(skipna=True)
    thresholds   = np.full(12, np.nan)
    rdi_vals     = pd.Series(np.nan, index=monthly_s.index)

    for m in range(1, 13):
        mask            = monthly_s.index.month == m
        q25             = monthly_s[mask].quantile(0.25)
        thresholds[m-1] = q25
        rdi_vals[mask]  = (monthly_s[mask] - q25) / overall_mean

    return rdi_vals, thresholds


def aggregate_huc4(storage_mat: pd.DataFrame) -> pd.Series:
    """
    Aggregate-first approach, intersection months only:
        1. Keep only months where ALL reservoirs have valid data
        2. Sum storage across reservoirs: S_HUC4_t = sum_i(S_i_t)
        3. Compute RDI on the combined series
    """
    valid_mask        = storage_mat.notna().all(axis=1)
    storage_intersect = storage_mat[valid_mask]
    combined          = storage_intersect.sum(axis=1).replace(0, np.nan)
    rdi_huc4, _       = compute_rdi(combined)
    return rdi_huc4


# ── Main loop over all valid HUC4s ────────────────────────────────────────────
print("\n" + "=" * 60)
print("Computing RDI per HUC4 ...")
print("=" * 60)

huc4_rows  = []   # one row per (HUC4, month)
resv_rows  = []   # one row per (reservoir, month)
skipped    = []   # HUC4s with no HydroShare data

for _, huc_row in valid_huc4.iterrows():
    huc4 = huc_row["HUC4"]
    name = huc_row["HUC4_Name"]

    # Reservoirs in this HUC4 that have HydroShare files
    resv = all_joined[
        (all_joined["huc4"] == huc4) & (all_joined["has_hs"])
    ].reset_index(drop=True)

    if resv.empty:
        skipped.append(huc4)
        print(f"  {huc4} {name:35s}  SKIP (no HydroShare data)")
        continue

    print(f"  {huc4} {name:35s}  {len(resv)} reservoirs", end="", flush=True)

    # Per-reservoir: load -> monthly avg
    storage_cols = {}

    for _, r in resv.iterrows():
        gid = int(r["GRanD_ID"])
        cap = r["max_historical_Storage_af"]
        try:
            daily   = load_hydroshare(gid, cap)
            monthly = daily.resample("MS").mean()

            storage_cols[gid] = monthly

            # Save per-reservoir rows (full individual range)
            for dt, s_val in zip(monthly.index, monthly.values):
                resv_rows.append({
                    "huc4": huc4, "huc4_name": name,
                    "GRanD_ID": gid, "dam_name": r["dam name"],
                    "capacity_af": cap,
                    "date": dt,
                    "Year": dt.year, "Month": dt.month,
                    "avg_storage_af": s_val,
                })
        except Exception as e:
            print(f"\n    WARNING: {r['dam name']} (ID {gid}) failed: {e}")

    if not storage_cols:
        skipped.append(huc4)
        print("  -> all reservoirs failed, skip")
        continue

    # HUC4 aggregation: intersection months only, then compute RDI
    storage_mat = pd.DataFrame(storage_cols)
    rdi_huc4    = aggregate_huc4(storage_mat)
    n_resv      = len(storage_cols)

    if rdi_huc4.empty:
        skipped.append(huc4)
        print("  -> no overlapping months across reservoirs, skip")
        continue

    for dt, rdi_val in zip(rdi_huc4.index, rdi_huc4.values):
        huc4_rows.append({
            "huc4": huc4, "huc4_name": name,
            "date": dt,
            "Year": dt.year, "Month": dt.month,
            "RDI_HUC4": rdi_val,
            "n_reservoirs": n_resv,
        })

    valid_n = sum(~np.isnan(rdi_huc4.values))
    if valid_n == 0:
        print(f"  valid=0mo  RDI=[n/a]")
    else:
        rng = (float(np.nanmin(rdi_huc4)), float(np.nanmax(rdi_huc4)))
        print(f"  valid={valid_n}mo  RDI=[{rng[0]:.2f},{rng[1]:.2f}]")

# ── Save outputs ──────────────────────────────────────────────────────────────
print("\n" + "=" * 60)
print("Saving outputs ...")
print("=" * 60)

huc4_df = pd.DataFrame(huc4_rows)
resv_df = pd.DataFrame(resv_rows)

p1 = os.path.join(OUTPUT_DIR, "all_HUC4_RDI.csv")
huc4_df.to_csv(p1, index=False, encoding="utf-8-sig")
print(f"  HUC4 RDI table     : {p1}  ({len(huc4_df)} rows)")

p2 = os.path.join(OUTPUT_DIR, "all_reservoir_monthly.csv")
resv_df.to_csv(p2, index=False, encoding="utf-8-sig")
print(f"  Per-reservoir table: {p2}  ({len(resv_df)} rows)")

print(f"\nSummary:")
print(f"  HUC4s computed : {huc4_df['huc4'].nunique()}")
print(f"  HUC4s skipped  : {len(skipped)} -> {skipped}")
print(f"  Reservoirs used: {resv_df['GRanD_ID'].nunique()}")

print("\n--- Sample output (HUC4 0513, first 6 months) ---")
sample = huc4_df[huc4_df["huc4"] == "0513"].head(6)
if not sample.empty:
    print(sample[["date","Year","Month","RDI_HUC4","n_reservoirs"]].to_string(index=False))

print("\n" + "=" * 60)
print("Done.")
print("=" * 60)
