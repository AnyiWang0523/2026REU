"""
select_ssi_gages.py  -  Identify the real USGS gage to use for SSI per HUC4
============================================================================
Method :
    - Use the gage closest to the point where the HUC4's main stem enters
      the HUC4 from upstream ("inlet"), NOT the outlet gage. Rationale:
      an outlet gage sits downstream of any reservoirs in the HUC4, so its
      flow already reflects reservoir release decisions -- exactly what
      WSDI is meant to help a manager decide whether to change. Using the
      inlet keeps SSI as an upstream/natural-supply signal, independent of
      in-HUC4 reservoir operations.
    - Headwater HUC4s (no upstream inflow -- the river originates inside
      the HUC4) are skipped entirely: their local tributaries are rainfall-
      driven and already captured by SPEI, so no SSI is computed for them.
    - HUC4s with multiple inlets (multiple upstream tributaries) keep all
      inlet gages; their flows are aggregated later (aggregate-first, same
      pattern as compute_rdi.py) before computing SSI.

Network topology source:
    data/Major_Streams/rs16my07.shp  -- EPA Reach File 1 (RF1)
    Each reach carries HUC/SEG (its own HUC8 + segment id) and
    DSHUC/DSSEG (the HUC8 + segment id it flows into). A reach whose own
    HUC4 differs from its downstream reach's HUC4 is a boundary-crossing
    reach; the reach it flows into (inside the target HUC4) is the inlet
    reach for that HUC4.

Gage source:
    data/CSV/USGS_Gages.csv -- 16,903 USGS gage locations (lat/lon).
    NOTE: this file only has a long-term mean discharge, not a monthly
    time series. Fetching the actual time series for the selected gages
    is a separate, later step.

Output:
    output/ssi/huc4_inlet_gages.csv  -- one row per (HUC4, inlet reach),
        with the nearest real USGS gage and distance in km.
    output/ssi/huc4_headwater_skipped.csv -- HUC4s with no upstream inlet.
"""

import os
import numpy as np
import pandas as pd
import geopandas as gpd

BASE_DIR   = r"C:\Users\anyiw\OneDrive - University of Illinois - Urbana\GIS\2026REU\WSDI"
VALID_CSV  = os.path.join(BASE_DIR, "output", "huc4_valid_coverage50.csv")
REACH_SHP  = os.path.join(BASE_DIR, "data", "Major_Streams", "rs16my07.shp")
GAGES_CSV  = os.path.join(BASE_DIR, "data", "CSV", "USGS_Gages.csv")
OUTPUT_DIR = os.path.join(BASE_DIR, "output", "ssi")
os.makedirs(OUTPUT_DIR, exist_ok=True)

# Projected CRS for accurate distance calculations (CONUS Albers Equal Area)
METRIC_CRS = "EPSG:5070"


def huc8_to_huc4(huc8) -> str:
    """First 4 digits of a zero-padded 8-digit HUC code."""
    return str(int(huc8)).zfill(8)[:4]


print("=" * 60)
print("Loading data ...")
print("=" * 60)

valid_huc4 = pd.read_csv(VALID_CSV, dtype={"HUC4": str})
valid_huc4["HUC4"] = valid_huc4["HUC4"].str.zfill(4)
target_huc4s = set(valid_huc4["HUC4"])
print(f"  Target HUC4s (from huc4_valid_coverage50.csv): {len(target_huc4s)}")

reaches = gpd.read_file(REACH_SHP)
reaches["HUC4_own"] = reaches["HUC"].apply(huc8_to_huc4)
reaches["HUC4_ds"]  = reaches["DSHUC"].apply(lambda d: huc8_to_huc4(d) if d and int(d) != 0 else None)
print(f"  Reaches loaded: {len(reaches)}")

gages = pd.read_csv(GAGES_CSV, dtype={"site_no": str})
gages_gdf = gpd.GeoDataFrame(
    gages, crs="EPSG:4326",
    geometry=gpd.points_from_xy(gages["dec_long_va"], gages["dec_lat_va"]),
).to_crs(METRIC_CRS)
print(f"  USGS gages loaded: {len(gages_gdf)}")

# Fast lookup: (HUC, SEG) -> reach row index, to find the receiving reach
reach_key_to_idx = {(row.HUC, row.SEG): i for i, row in enumerate(reaches.itertuples())}

# ── Find boundary-crossing (inlet-defining) reaches for each target HUC4 ──
print("\n" + "=" * 60)
print("Finding inlet reaches per HUC4 ...")
print("=" * 60)

reaches_metric = reaches.to_crs(METRIC_CRS)

inlet_rows    = []
headwater_rows = []

for huc4 in sorted(target_huc4s):
    name = valid_huc4.loc[valid_huc4["HUC4"] == huc4, "HUC4_Name"].iloc[0]

    # Reaches upstream of (or crossing into) this HUC4: their own HUC4 differs
    # from this target, but their downstream HUC4 equals this target.
    crossing = reaches[(reaches["HUC4_own"] != huc4) & (reaches["HUC4_ds"] == huc4)]

    inlet_reach_idxs = set()
    for row in crossing.itertuples():
        key = (row.DSHUC, row.DSSEG)
        if key in reach_key_to_idx:
            inlet_reach_idxs.add(reach_key_to_idx[key])

    if not inlet_reach_idxs:
        headwater_rows.append({"HUC4": huc4, "HUC4_Name": name})
        print(f"  {huc4} {name:35s}  HEADWATER (no upstream inlet) -> skip SSI")
        continue

    print(f"  {huc4} {name:35s}  {len(inlet_reach_idxs)} inlet reach(es)")

    for idx in inlet_reach_idxs:
        inlet_geom = reaches_metric.geometry.iloc[idx]
        dists = gages_gdf.geometry.distance(inlet_geom)
        nearest_i = dists.values.argmin()
        nearest_gage = gages_gdf.iloc[nearest_i]
        dist_km = dists.values[nearest_i] / 1000.0

        inlet_rows.append({
            "HUC4": huc4,
            "HUC4_Name": name,
            "inlet_HUC": reaches.iloc[idx]["HUC"],
            "inlet_SEG": reaches.iloc[idx]["SEG"],
            "site_no": nearest_gage["site_no"],
            "station_nm": nearest_gage["station_nm"],
            "gage_lat": nearest_gage["dec_lat_va"],
            "gage_lon": nearest_gage["dec_long_va"],
            "distance_km": round(dist_km, 2),
        })

# ── Save outputs ──────────────────────────────────────────────────────────
inlet_df    = pd.DataFrame(inlet_rows)
headwater_df = pd.DataFrame(headwater_rows)

p1 = os.path.join(OUTPUT_DIR, "huc4_inlet_gages.csv")
inlet_df.to_csv(p1, index=False, encoding="utf-8-sig")

p2 = os.path.join(OUTPUT_DIR, "huc4_headwater_skipped.csv")
headwater_df.to_csv(p2, index=False, encoding="utf-8-sig")

print("\n" + "=" * 60)
print("Summary")
print("=" * 60)
print(f"  HUC4s with inlet gage(s) : {inlet_df['HUC4'].nunique()} ({len(inlet_df)} inlet reaches total)")
print(f"  HUC4s skipped (headwater): {len(headwater_df)}")
print(f"  -> {p1}")
print(f"  -> {p2}")
