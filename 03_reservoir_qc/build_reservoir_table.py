# ============================================================
# Reservoir Coverage Table Builder
# Columns:
#   1. HUC4 Name
#   2. GDROM reservoirs lat/lon list
#   3. GDROM total max storage (acre-feet) per HUC4
#   4. NID (>1km2) reservoirs lat/lon list (deduplicated) + total max storage
#   5a. GDROM total / NID total (%) — paper Criteria 2 metric
#   5b. GDROM total / NID name-matched (%) — diagnostic only
# ============================================================

import geopandas as gpd
import pandas as pd
import numpy as np
import os
from difflib import SequenceMatcher
from math import radians, sin, cos, sqrt, atan2

BASE_DIR   = r"C:\Users\anyiw\OneDrive - University of Illinois - Urbana\GIS\2026REU\WSDI"
HUC4_SHP   = os.path.join(BASE_DIR, "data", "HUC4", "WBDHU4.shp")
NID_SHP    = os.path.join(BASE_DIR, "data", "NID_All_US_over1km2", "NID_All_USA_Over_1km2.shp")
GDROM_CSV  = os.path.join(BASE_DIR, "data", "GDROM_Resv.csv")
OUTPUT     = os.path.join(BASE_DIR, "output", "reservoir_coverage_table.csv")
OUTPUT_SEL = os.path.join(BASE_DIR, "output", "reservoir_coverage_selected.csv")
os.makedirs(os.path.join(BASE_DIR, "output"), exist_ok=True)

# ============================================================
# STEP 1: Load data
# ============================================================
print("=" * 55)
print("STEP 1: Loading data...")
print("=" * 55)

huc4_gdf = gpd.read_file(HUC4_SHP)
huc4_gdf["huc4"] = huc4_gdf["huc4"].astype(str).str.zfill(4)
print(f"  HUC4 polygons: {len(huc4_gdf)}")

gdrom_df = pd.read_csv(GDROM_CSV).rename(columns={
    "lattitude":  "latitude",
    "longtitude": "longitude",
})
gdrom_df = gdrom_df.dropna(subset=["latitude", "longitude"])
print(f"  GDROM reservoirs: {len(gdrom_df)}")

nid_gdf = gpd.read_file(NID_SHP)
nid_gdf["maxStorage_num"] = pd.to_numeric(nid_gdf["maxStorage"], errors="coerce")
nid_gdf["latitude_num"]   = pd.to_numeric(nid_gdf["latitude"],   errors="coerce")
nid_gdf["longitude_num"]  = pd.to_numeric(nid_gdf["longitude"],  errors="coerce")
nid_gdf = nid_gdf.dropna(subset=["latitude_num", "longitude_num", "maxStorage_num"])
nid_gdf = nid_gdf[nid_gdf["maxStorage_num"] > 0].copy()
print(f"  NID >1km2 reservoirs (valid): {len(nid_gdf)}")

# ============================================================
# STEP 2: Spatial join GDROM -> HUC4
# ============================================================
print("\n" + "=" * 55)
print("STEP 2: Spatial join GDROM -> HUC4...")
print("=" * 55)

gdrom_gdf = gpd.GeoDataFrame(
    gdrom_df,
    geometry=gpd.points_from_xy(gdrom_df["longitude"], gdrom_df["latitude"]),
    crs="EPSG:4326"
).to_crs(huc4_gdf.crs)

gdrom_joined = gpd.sjoin(
    gdrom_gdf[["GRanD_ID", "dam name", "latitude", "longitude",
               "max_historical_Storage_af", "geometry"]],
    huc4_gdf[["huc4", "name", "geometry"]],
    how="left",
    predicate="within"
)
gdrom_joined = gdrom_joined.rename(columns={"name": "HUC4_Name"})
matched = gdrom_joined["huc4"].notna().sum()
print(f"  GDROM matched to HUC4: {matched} / {len(gdrom_joined)}")

# ============================================================
# STEP 3: Spatial join NID -> HUC4
# ============================================================
print("\n" + "=" * 55)
print("STEP 3: Spatial join NID -> HUC4...")
print("=" * 55)

nid_proj = nid_gdf.to_crs(huc4_gdf.crs)

nid_joined = gpd.sjoin(
    nid_proj[["name", "latitude_num", "longitude_num", "maxStorage_num", "geometry"]],
    huc4_gdf[["huc4", "name", "geometry"]],
    how="left",
    predicate="within"
).rename(columns={"name_left": "dam_name", "name_right": "HUC4_Name"})

matched_nid = nid_joined["huc4"].notna().sum()
print(f"  NID matched to HUC4: {matched_nid} / {len(nid_joined)}")

# ============================================================
# STEP 4: Deduplicate NID within each HUC4
#         Rule: same maxStorage AND lat/lon within 0.01 deg -> same reservoir
# ============================================================
print("\n" + "=" * 55)
print("STEP 4: Deduplicating NID reservoirs...")
print("=" * 55)

def dedup_nid_group(group):
    """
    Within a HUC4, deduplicate reservoirs that are 'the same':
    same maxStorage AND lat/lon within 0.01 degree of each other.
    Returns deduplicated DataFrame.
    """
    if group.empty:
        return group
    rows = group.reset_index(drop=True)
    keep = [True] * len(rows)
    for i in range(len(rows)):
        if not keep[i]:
            continue
        for j in range(i + 1, len(rows)):
            if not keep[j]:
                continue
            same_storage = abs(rows.loc[i, "maxStorage_num"] - rows.loc[j, "maxStorage_num"]) < 1e-3
            dlat = abs(rows.loc[i, "latitude_num"]  - rows.loc[j, "latitude_num"])
            dlon = abs(rows.loc[i, "longitude_num"] - rows.loc[j, "longitude_num"])
            if same_storage and dlat < 0.01 and dlon < 0.01:
                keep[j] = False
    return rows[keep]

nid_valid = nid_joined.dropna(subset=["huc4"]).copy()
nid_dedup_parts = []
for _, grp in nid_valid.groupby("huc4"):
    nid_dedup_parts.append(dedup_nid_group(grp))
nid_dedup = pd.concat(nid_dedup_parts).reset_index(drop=True)
print(f"  NID before dedup: {len(nid_valid)}")
print(f"  NID after dedup:  {len(nid_dedup)}")

# ============================================================
# STEP 5: Build per-HUC4 summary
# ============================================================
print("\n" + "=" * 55)
print("STEP 5: Building HUC4 summary table...")
print("=" * 55)

def fmt_coords(sub, lat_col, lon_col):
    """Format list of (lat, lon) as readable string."""
    pairs = [f"({row[lat_col]:.4f}, {row[lon_col]:.4f})"
             for _, row in sub.iterrows()]
    return "; ".join(pairs) if pairs else ""

def haversine_km(lat1, lon1, lat2, lon2):
    """Great-circle distance in km between two (lat, lon) points."""
    R = 6371.0
    phi1, phi2 = radians(lat1), radians(lat2)
    dphi = radians(lat2 - lat1)
    dlam = radians(lon2 - lon1)
    a = sin(dphi / 2)**2 + cos(phi1) * cos(phi2) * sin(dlam / 2)**2
    return R * 2 * atan2(sqrt(a), sqrt(1 - a))

def find_best_nid_match(gdrom_name, gdrom_lat, gdrom_lon, gdrom_storage_af,
                        nid_grp,
                        name_threshold=0.6, max_dist_km=80.0,
                        nearby_km=5.0, stor_ratio_min=0.1, stor_ratio_max=10.0):
    """
    Match a GDROM dam to a NID dam by two OR conditions:
      1. Name: name_similarity >= 0.6 within 80 km  (original)
      2. Proximity + storage: dist < 5 km AND storage ratio in [0.1, 10]
         (catches co-located dams with different name spellings)
    Returns (best_row, name_score, dist_km) or (None, 0, None).
    """
    gdrom_clean = str(gdrom_name).lower().strip()

    # Collect all candidates within max_dist_km
    candidates = []
    for _, row in nid_grp.iterrows():
        dist_km = haversine_km(gdrom_lat, gdrom_lon,
                               row["latitude_num"], row["longitude_num"])
        if dist_km > max_dist_km:
            continue
        name_score = SequenceMatcher(None, gdrom_clean,
                                     str(row["dam_name"]).lower().strip()).ratio()
        prox_score = max(0.0, 1.0 - dist_km / max_dist_km)
        combined   = 0.6 * name_score + 0.4 * prox_score
        candidates.append((combined, name_score, dist_km, row))

    if not candidates:
        return None, 0.0, None

    candidates.sort(key=lambda x: x[0], reverse=True)
    _, best_name_score, best_dist, best_row = candidates[0]

    # Condition 1: high name similarity
    if best_name_score >= name_threshold:
        return best_row, best_name_score, best_dist

    # Condition 2: very close location + reasonable storage ratio
    for _, name_score, dist_km, row in candidates:
        if dist_km > nearby_km:
            continue
        n_stor = row["maxStorage_num"]
        if n_stor > 0 and gdrom_storage_af > 0:
            ratio = gdrom_storage_af / n_stor
            if stor_ratio_min <= ratio <= stor_ratio_max:
                return row, name_score, dist_km  # proximity+storage match

    return None, 0.0, None

gdrom_valid = gdrom_joined.dropna(subset=["huc4"])

gdrom_summary = []
suspicious_pairs = []   # dams flagged for review (no NID match, or extreme storage ratio)

for huc4_id, grp in gdrom_valid.groupby("huc4"):
    huc4_name = grp["HUC4_Name"].iloc[0]
    coords    = fmt_coords(grp, "latitude", "longitude")
    names     = "; ".join(grp["dam name"].fillna("").astype(str).tolist())
    total_str = grp["max_historical_Storage_af"].sum()

    nid_grp = nid_dedup[nid_dedup["huc4"] == huc4_id]
    matched_nid_storage = 0.0
    for _, gdrom_row in grp.iterrows():
        match, name_sc, dist_km = find_best_nid_match(
            gdrom_row["dam name"],
            gdrom_row["latitude"], gdrom_row["longitude"],
            gdrom_row["max_historical_Storage_af"],
            nid_grp
        )
        if match is None:
            # No NID counterpart found — flag for review
            suspicious_pairs.append({
                "flag_type":  "no_nid_match",
                "huc4":       huc4_id,
                "GDROM_name": gdrom_row["dam name"],
                "NID_name":   "",
                "name_score": np.nan,
                "dist_km":    np.nan,
                "GDROM_af":   round(gdrom_row["max_historical_Storage_af"], 0),
                "NID_af":     np.nan,
                "ratio_G_N":  np.nan,
                "note":       "No NID match: name_sim < 0.6 within 80 km, and no co-located (< 5 km) NID dam with similar storage",
            })
        else:
            matched_nid_storage += match["maxStorage_num"]
            # Flag pairs where storage differs by >5x (potential definition mismatch)
            g_stor = gdrom_row["max_historical_Storage_af"]
            n_stor = match["maxStorage_num"]
            if g_stor > 0 and n_stor > 0:
                ratio = g_stor / n_stor
                if ratio > 5 or ratio < 0.2:
                    suspicious_pairs.append({
                        "flag_type":  "storage_ratio",
                        "huc4":       huc4_id,
                        "GDROM_name": gdrom_row["dam name"],
                        "NID_name":   match["dam_name"],
                        "name_score": round(name_sc, 3),
                        "dist_km":    round(dist_km, 2),
                        "GDROM_af":   round(g_stor, 0),
                        "NID_af":     round(n_stor, 0),
                        "ratio_G_N":  round(ratio, 2),
                        "note": ("GDROM>>NID: check if NID uses normal storage"
                                 if ratio > 5 else
                                 "NID>>GDROM: NID design capacity >> observed max (e.g. flood-control dam)"),
                    })

    gdrom_summary.append({
        "huc4":                      huc4_id,
        "HUC4_Name":                 huc4_name,
        "GDROM_Coords":              coords,
        "GDROM_Names":               names,
        "GDROM_Count":               len(grp),
        "GDROM_TotalStorage_af":     round(total_str, 2),
        "NID_NameMatched_Storage_af": round(matched_nid_storage, 2) if matched_nid_storage > 0 else np.nan
    })

gdrom_sum_df = pd.DataFrame(gdrom_summary)

# ---- NID per HUC4 ----
nid_summary = []
for huc4_id, grp in nid_dedup.groupby("huc4"):
    huc4_name = grp["HUC4_Name"].iloc[0]
    coords    = fmt_coords(grp, "latitude_num", "longitude_num")
    total_str = grp["maxStorage_num"].sum()
    nid_summary.append({
        "huc4":                huc4_id,
        "NID_Coords":          coords,
        "NID_Count":           len(grp),
        "NID_TotalStorage_af": round(total_str, 2)
    })

nid_sum_df = pd.DataFrame(nid_summary)

# ---- Merge with HUC4 name list ----
huc4_names = huc4_gdf[["huc4", "name"]].rename(columns={"name": "HUC4_Name"})

result = (
    huc4_names
    .merge(gdrom_sum_df[["huc4", "GDROM_Coords", "GDROM_Names", "GDROM_Count",
                          "GDROM_TotalStorage_af", "NID_NameMatched_Storage_af"]],
           on="huc4", how="left")
    .merge(nid_sum_df[["huc4", "NID_Coords", "NID_Count", "NID_TotalStorage_af"]],
           on="huc4", how="left")
)

# ---- Column 5: NID matched / NID total — both numerator and denominator from NID ----
# Answers: what fraction of total NID capacity in this HUC4 was matched to a GDROM reservoir?
raw_pct = np.where(
    (result["NID_TotalStorage_af"] > 0) & (result["NID_NameMatched_Storage_af"].notna()),
    (result["NID_NameMatched_Storage_af"] / result["NID_TotalStorage_af"] * 100),
    np.nan
)
result["Coverage_Pct"] = np.where(
    ~np.isnan(raw_pct.astype(float)),
    raw_pct.astype(float).round(2),
    np.nan
)

# ---- Diagnostic: GDROM / NID total (different data sources, kept for reference) ----
raw_pct_gdrom = np.where(
    (result["NID_TotalStorage_af"] > 0) & (result["GDROM_TotalStorage_af"].notna()),
    (result["GDROM_TotalStorage_af"] / result["NID_TotalStorage_af"] * 100),
    np.nan
)
result["Coverage_Pct_GDROM"] = np.where(
    ~np.isnan(raw_pct_gdrom.astype(float)),
    raw_pct_gdrom.astype(float).round(2),
    np.nan
)

# ---- Fill NaN counts with 0 ----
result["GDROM_Count"] = result["GDROM_Count"].fillna(0).astype(int)
result["NID_Count"]   = result["NID_Count"].fillna(0).astype(int)

# ---- Criteria 2 pass flag: 50% <= Col5 <= 100% ----
result["Criteria2_Pass"] = (
    result["Coverage_Pct"].notna() &
    result["Coverage_Pct"].between(50, 100)
)

# ---- Final column order ----
result = result[[
    "huc4",
    "HUC4_Name",
    "GDROM_Coords",
    "GDROM_Names",
    "GDROM_Count",
    "GDROM_TotalStorage_af",
    "NID_Coords",
    "NID_Count",
    "NID_TotalStorage_af",
    "NID_NameMatched_Storage_af",
    "Coverage_Pct",
    "Coverage_Pct_GDROM",
    "Criteria2_Pass",
]]

result.columns = [
    "HUC4",
    "HUC4_Name",
    "Col2_GDROM_Reservoir_LatLon",
    "Col2b_GDROM_Reservoir_Names",
    "Col2c_GDROM_Count",
    "Col3_GDROM_Total_MaxStorage_af",
    "Col4_NID_Reservoir_LatLon_Deduped",
    "Col4b_NID_Count",
    "Col4c_NID_Total_MaxStorage_af",
    "Col4d_NID_NameMatched_MaxStorage_af",
    "Col5_Coverage_Pct_NIDMatched_vs_NIDTotal",   # Criteria 2: NID/NID, use for 50% filter
    "Col5b_Coverage_Pct_GDROM_vs_NIDTotal",        # diagnostic: GDROM/NID
    "Criteria2_Pass",
]

# ============================================================
# STEP 5b: Flag HUC4s where Col5b > 100%
# ============================================================
print("\n" + "=" * 55)
print("STEP 5b: Flagging HUC4s with Col5b > 100%...")
print("=" * 55)

over100 = result[
    result["Col5_Coverage_Pct_NIDMatched_vs_NIDTotal"].notna() &
    (result["Col5_Coverage_Pct_NIDMatched_vs_NIDTotal"] > 100)
]
for _, row in over100.iterrows():
    suspicious_pairs.append({
        "flag_type":  "coverage_over100",
        "huc4":       row["HUC4"],
        "GDROM_name": row["Col2b_GDROM_Reservoir_Names"],
        "NID_name":   "",
        "name_score": np.nan,
        "dist_km":    np.nan,
        "GDROM_af":   row["Col3_GDROM_Total_MaxStorage_af"],
        "NID_af":     row["Col4d_NID_NameMatched_MaxStorage_af"],
        "ratio_G_N":  np.nan,
        "note":       f"HUC4 Col5 = {row['Col5_Coverage_Pct_NIDMatched_vs_NIDTotal']:.1f}%: NID matched exceeds NID total (dedup issue)",
    })
print(f"  HUC4s with Col5b > 100%: {len(over100)}")

# ============================================================
# STEP 6: Save
# ============================================================
print("\n" + "=" * 55)
print("STEP 6: Saving output...")
print("=" * 55)

result.to_csv(OUTPUT, index=False, encoding="utf-8-sig")
print(f"  Saved full table: {OUTPUT}")

# Save all flagged dams / HUC4s for manual inspection (before selected, so Excel lock can't block it)
if suspicious_pairs:
    susp_path = os.path.join(BASE_DIR, "output", "reservoir_storage_flags.csv")
    flags_df = pd.DataFrame(suspicious_pairs).sort_values(
        ["flag_type", "ratio_G_N"], ascending=[True, False], na_position="last"
    )
    flags_df.to_csv(susp_path, index=False, encoding="utf-8-sig")
    no_match_n     = (flags_df["flag_type"] == "no_nid_match").sum()
    ratio_flag_n   = (flags_df["flag_type"] == "storage_ratio").sum()
    over100_flag_n = (flags_df["flag_type"] == "coverage_over100").sum()
    print(f"  Flagged {len(suspicious_pairs)} entries -> {susp_path}")
    print(f"    no_nid_match    : {no_match_n}  (dam level)")
    print(f"    storage_ratio   : {ratio_flag_n}  (dam level, ratio >5x or <0.2x)")
    print(f"    coverage_over100: {over100_flag_n}  (HUC4 level, Col5b > 100%)")

print(f"  Total HUC4 rows: {len(result)}")
print(f"  HUC4 with GDROM data:        {(result['Col3_GDROM_Total_MaxStorage_af'] > 0).sum()}")
print(f"  HUC4 with NID data:          {(result['Col4c_NID_Total_MaxStorage_af'] > 0).sum()}")
print(f"  HUC4 with name-matched NID:  {result['Col4d_NID_NameMatched_MaxStorage_af'].notna().sum()}")
print(f"  HUC4 with NID/NID coverage : {result['Col5_Coverage_Pct_NIDMatched_vs_NIDTotal'].notna().sum()}")
print(f"  HUC4 with GDROM/NID coverage:{result['Col5b_Coverage_Pct_GDROM_vs_NIDTotal'].notna().sum()}")
print(f"  HUC4 passing Criteria 2:     {result['Criteria2_Pass'].sum()}")

print("\n--- Preview (HUC4 with GDROM data) ---")
preview = result[result["Col3_GDROM_Total_MaxStorage_af"].notna()].head(5)
for _, row in preview.iterrows():
    print(f"\n  HUC4: {row['HUC4']} | {row['HUC4_Name']}")
    print(f"    GDROM storage:        {row['Col3_GDROM_Total_MaxStorage_af']:,.0f} af")
    print(f"    NID total storage:    {row['Col4c_NID_Total_MaxStorage_af']:,.0f} af" if pd.notna(row['Col4c_NID_Total_MaxStorage_af']) else "    NID total storage: N/A")
    print(f"    NID matched storage:  {row['Col4d_NID_NameMatched_MaxStorage_af']:,.0f} af" if pd.notna(row['Col4d_NID_NameMatched_MaxStorage_af']) else "    NID matched storage: N/A")
    print(f"    Coverage NID/NID:    {row['Col5_Coverage_Pct_NIDMatched_vs_NIDTotal']}%" if pd.notna(row['Col5_Coverage_Pct_NIDMatched_vs_NIDTotal']) else "    Coverage NID/NID: N/A")
    print(f"    Coverage GDROM/NID: {row['Col5b_Coverage_Pct_GDROM_vs_NIDTotal']}%" if pd.notna(row['Col5b_Coverage_Pct_GDROM_vs_NIDTotal']) else "    Coverage GDROM/NID: N/A")
    print(f"    Criteria2 Pass:       {row['Criteria2_Pass']}")

# Save HUC4s that pass Criteria 2 (last — OK if Excel has it open, non-critical)
selected = result[result["Criteria2_Pass"]]
try:
    selected.to_csv(OUTPUT_SEL, index=False, encoding="utf-8-sig")
    print(f"  Saved Criteria2-selected: {OUTPUT_SEL}  ({len(selected)} HUC4s)")
except PermissionError:
    print(f"  WARNING: could not write {OUTPUT_SEL} (file open in Excel?). Non-critical — filter_valid_huc4.py reads the full table.")

print("\n" + "=" * 55)
print("Done!")
print("=" * 55)
