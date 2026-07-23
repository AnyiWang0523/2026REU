"""
compute_ssi.py  -  SSI-3 (Standardized Streamflow Index) for selected inlet gages
==================================================================================
Method (Shukla & Wood, 2008): identical procedure to SPI, substituting monthly
mean streamflow for precipitation -- reuses compute_spi() from
src/compute_indices.py (accumulate -> gamma fit -> probit transform), applied
to discharge instead of rainfall.

Gage source:
    output/ssi/huc4_inlet_gages.csv  -- inlet gage(s) per HUC4, from
    05_streamflow/select_ssi_gages.py.

Data source:
    USGS NWIS daily discharge (parameterCd 00060), fetched via
    dataretrieval.nwis.get_dv for the fixed 2000-01-01..2020-12-31 window --
    the same calibration window used for SPI/SPEI in
    04_drought_analysis/drought_index_analysis.py (CALIB_START/END = 2000, 2020).

    Not every selected inlet gage has data in this window (checked against the
    same NWIS backend that https://apps.usgs.gov/nwismapper/ queries -- some
    gages only cover e.g. 2012-2026, others only 1974-1993). Those gages are
    skipped and written to gages_no_2000_2020_data.csv.

HUC4 aggregation -- aggregate-first (same pattern as 03_reservoir_qc/compute_rdi.py):
    when a HUC4 has multiple inlet gages, monthly flows are summed across
    gages over the intersection of months where ALL of that HUC4's inlet
    gages have data, then SSI-3 is computed on the combined series.

Outputs:
    output/ssi/HUC4_SSI3.csv               -- monthly SSI-3 per valid HUC4
    output/ssi/gages_no_2000_2020_data.csv -- gages with no discharge data
                                               in 2000-2020 (written out, per
                                               request, instead of silently
                                               dropped)
"""

import os
import sys
import time
import warnings
import numpy as np
import pandas as pd
import dataretrieval.nwis as nwis

warnings.filterwarnings("ignore")

BASE_DIR = r"C:\Users\anyiw\OneDrive - University of Illinois - Urbana\GIS\2026REU\WSDI"
sys.path.insert(0, BASE_DIR)
from src.compute_indices import compute_spi

GAGES_CSV  = os.path.join(BASE_DIR, "output", "ssi", "huc4_inlet_gages.csv")
OUTPUT_DIR = os.path.join(BASE_DIR, "output", "ssi")
os.makedirs(OUTPUT_DIR, exist_ok=True)

START, END             = "2000-01-01", "2020-12-31"
CALIB_START, CALIB_END = 2000, 2020
FULL_MONTHS            = pd.date_range("2000-01-01", "2020-12-01", freq="MS")

# Minimum valid months required before attempting a per-calendar-month gamma
# fit for a HUC4's combined flow series (fewer than this gives too few
# samples per calendar month to trust the fit). Judgment call, not a fixed
# rule -- adjust if too many/too few HUC4s end up excluded.
MIN_VALID_MONTHS = 24


def fetch_monthly_flow(site_no: str):
    """
    Fetch daily discharge for `site_no` over 2000-2020, return monthly mean
    (cfs) reindexed to the full 2000-01..2020-12 month range (NaN where
    missing). Returns None if the site has no data at all in this window.
    """
    try:
        df, _ = nwis.get_dv(sites=site_no, parameterCd="00060", start=START, end=END)
    except Exception as e:
        print(f"    WARNING: {site_no} query failed: {e}")
        return None

    if df.empty:
        return None

    flow_col = [c for c in df.columns if "00060" in c and "Mean" in c]
    if not flow_col:
        flow_col = [c for c in df.columns if "00060" in c]
    if not flow_col:
        return None

    s = pd.to_numeric(df[flow_col[0]], errors="coerce")
    s.index = pd.to_datetime(df.index).tz_localize(None)   # NWIS returns tz-aware (UTC) index
    s = s.dropna()
    if s.empty:
        return None

    monthly = s.resample("MS").mean().reindex(FULL_MONTHS)
    return monthly


# ── Load selected inlet gages ─────────────────────────────────────────────
print("=" * 60)
print("Loading selected inlet gages ...")
print("=" * 60)

gages_df = pd.read_csv(GAGES_CSV, dtype={"HUC4": str, "site_no": str})
gages_df["HUC4"] = gages_df["HUC4"].str.zfill(4)
unique_sites = gages_df[["site_no", "station_nm"]].drop_duplicates()
print(f"  Inlet gage rows : {len(gages_df)}  ({len(unique_sites)} unique sites)")

# ── Fetch each unique gage once, check 2000-2020 coverage ─────────────────
print("\n" + "=" * 60)
print("Fetching 2000-2020 discharge per gage ...")
print("=" * 60)

flow_by_site = {}
no_data_rows = []

for _, row in unique_sites.iterrows():
    site, name = row["site_no"], row["station_nm"]
    monthly = fetch_monthly_flow(site)
    if monthly is None or monthly.notna().sum() == 0:
        print(f"  {site}  {name:45s}  NO DATA in 2000-2020")
        no_data_rows.append({"site_no": site, "station_nm": name})
    else:
        n_valid = int(monthly.notna().sum())
        print(f"  {site}  {name:45s}  {n_valid}/252 months")
        flow_by_site[site] = monthly
    time.sleep(0.3)

# ── Aggregate-first per HUC4 (sum inlet gages over shared valid months) ───
print("\n" + "=" * 60)
print("Aggregating to HUC4 and computing SSI-3 ...")
print("=" * 60)

huc4_flow_cols = {}
skipped_huc4   = []

for huc4, sub in gages_df.groupby("HUC4"):
    sites = [s for s in sub["site_no"] if s in flow_by_site]
    if not sites:
        skipped_huc4.append((huc4, "no usable gage data"))
        print(f"  {huc4}  SKIP (no usable gage data)")
        continue

    flow_mat   = pd.DataFrame({s: flow_by_site[s] for s in sites})
    valid_mask = flow_mat.notna().all(axis=1)
    combined   = flow_mat[valid_mask].sum(axis=1).reindex(FULL_MONTHS)

    n_valid = int(combined.notna().sum())
    if n_valid < MIN_VALID_MONTHS:
        skipped_huc4.append((huc4, f"only {n_valid} valid months"))
        print(f"  {huc4}  SKIP ({n_valid} valid months, below MIN_VALID_MONTHS={MIN_VALID_MONTHS})")
        continue

    huc4_flow_cols[huc4] = combined
    print(f"  {huc4}  {len(sites)} gage(s)  {n_valid}/252 valid months")

flow_df = pd.DataFrame(huc4_flow_cols, index=FULL_MONTHS)
ssi3_df = compute_spi(flow_df, scale=3, calib_start_year=CALIB_START, calib_end_year=CALIB_END)

# ── Save outputs ────────────────────────────────────────────────────────
print("\n" + "=" * 60)
print("Saving outputs ...")
print("=" * 60)

long_rows = []
for huc4 in ssi3_df.columns:
    for dt, val in ssi3_df[huc4].items():
        long_rows.append({"HUC4": huc4, "date": dt, "Year": dt.year, "Month": dt.month, "SSI3": val})

ssi_long = pd.DataFrame(long_rows, columns=["HUC4", "date", "Year", "Month", "SSI3"])
n_huc4_ok = ssi_long["HUC4"].nunique()
p1 = os.path.join(OUTPUT_DIR, "HUC4_SSI3.csv")
ssi_long.to_csv(p1, index=False, encoding="utf-8-sig")
print(f"  SSI-3 table          : {p1}  ({len(ssi_long)} rows, {n_huc4_ok} HUC4s)")

no_data_df = pd.DataFrame(no_data_rows, columns=["site_no", "station_nm"])
p2 = os.path.join(OUTPUT_DIR, "gages_no_2000_2020_data.csv")
no_data_df.to_csv(p2, index=False, encoding="utf-8-sig")
print(f"  Gages w/o 2000-2020  : {p2}  ({len(no_data_df)} gages)")

print(f"\nSummary:")
print(f"  HUC4s with SSI-3   : {n_huc4_ok}")
print(f"  HUC4s skipped      : {len(skipped_huc4)} -> {skipped_huc4}")
print(f"  Gages w/o data     : {len(no_data_df)} / {len(unique_sites)}")

print("\n" + "=" * 60)
print("Done.")
print("=" * 60)
