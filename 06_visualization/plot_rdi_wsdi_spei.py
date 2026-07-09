"""
Generate Figure: RDI / WSDI / SPEI-1 time series for one HUC4 (2000-2020).

Layout (3 rows, shared x-axis):
  Row 1 — SPEI-1      monthly    from output/drought_correlation/huc4_spei1.csv
  Row 2 — RDI         monthly    computed from USGS lake-elevation data (param 00054)
  Row 3 — WSDI        annual     from output/<WSDI_TYPE>/WSDI_*_HUC4_2000_2020.csv

HUC4 selection: auto-picks the first valid HUC4 that has data in all three
sources, in the order given by PREFERRED_HUC4S. Override with TARGET_HUC4.
"""

import os, sys, time
import numpy as np
# Force UTF-8 output so special characters don't crash on Windows cp1252 terminals
sys.stdout.reconfigure(encoding="utf-8", errors="replace")
import pandas as pd
import geopandas as gpd
import matplotlib.pyplot as plt
import matplotlib.dates as mdates
import matplotlib.patches as mpatches
import dataretrieval.nwis as nwis

# ── Configuration ─────────────────────────────────────────────────────────────
BASE_DIR  = r"C:\Users\anyiw\OneDrive - University of Illinois - Urbana\GIS\2026REU\WSDI"
WSDI_TYPE = "PS_GW"          # PS_GW | PS_SW | IR_GW | IR_SW

# Auto-select first HUC4 that has SPEI + WSDI + GDROM data.
# Set TARGET_HUC4 to a specific code (e.g. "1113") to override.
TARGET_HUC4     = None
PREFERRED_HUC4S = ["0513", "1113", "1403", "1014", "0512"]

# ── Paths ─────────────────────────────────────────────────────────────────────
VALID_CSV = os.path.join(BASE_DIR, "output", "huc4_valid_coverage50.csv")
GDROM_CSV = os.path.join(BASE_DIR, "data", "CSV", "GDROM_Resv.csv")
SPEI_CSV  = os.path.join(BASE_DIR, "output", "drought_correlation", "huc4_spei1.csv")
WSDI_CSV  = os.path.join(BASE_DIR, "output", WSDI_TYPE,
                          f"WSDI_{WSDI_TYPE}_HUC4_2000_2020.csv")
HUC4_SHP  = os.path.join(BASE_DIR, "data", "HUC4", "WBDHU4.shp")
OUT_DIR   = os.path.join(BASE_DIR, "output", "figures")
os.makedirs(OUT_DIR, exist_ok=True)

# ── Step 1: Choose target HUC4 ───────────────────────────────────────────────
print("=" * 55)
print("STEP 1: Selecting target HUC4...")
print("=" * 55)

valid  = pd.read_csv(VALID_CSV, dtype={"HUC4": str})
spei   = pd.read_csv(SPEI_CSV, index_col=0, parse_dates=True)
wsdi_df = pd.read_csv(WSDI_CSV, dtype={"HUC4": str})

valid_set = set(valid["HUC4"])
spei_set  = set(spei.columns)
wsdi_set  = set(wsdi_df["HUC4"].str.zfill(4))

candidates = [h for h in PREFERRED_HUC4S
              if h in valid_set and h in spei_set and h in wsdi_set]

if TARGET_HUC4 and TARGET_HUC4 in valid_set:
    huc4_id = TARGET_HUC4
elif candidates:
    huc4_id = candidates[0]
else:
    raise RuntimeError(
        f"No preferred HUC4 found in all three data sources.\n"
        f"SPEI covers: {sorted(spei_set)[:10]} ...\n"
        f"Valid HUC4s: {sorted(valid_set)[:10]} ..."
    )

huc4_name = valid.loc[valid["HUC4"] == huc4_id, "HUC4_Name"].values[0]
print(f"  Selected: {huc4_id} — {huc4_name}")

# ── Step 2: Find GDROM reservoirs in this HUC4 ───────────────────────────────
print("\n" + "=" * 55)
print("STEP 2: Loading GDROM reservoirs...")
print("=" * 55)

huc4_gdf = gpd.read_file(HUC4_SHP).to_crs("EPSG:4326")
huc4_gdf["huc4"] = huc4_gdf["huc4"].astype(str).str.zfill(4)
target_poly = huc4_gdf[huc4_gdf["huc4"] == huc4_id]

gdrom = pd.read_csv(GDROM_CSV).dropna(subset=["lattitude", "longtitude"])
gdrom["missing_ratio"] = pd.to_numeric(gdrom["missing_ratio"], errors="coerce")
gdrom_gdf = gpd.GeoDataFrame(
    gdrom,
    geometry=gpd.points_from_xy(gdrom["longtitude"], gdrom["lattitude"]),
    crs="EPSG:4326"
)
resv_in = gpd.sjoin(
    gdrom_gdf, target_poly[["huc4", "geometry"]], how="inner", predicate="within"
).reset_index(drop=True)

print(f"  Reservoirs in HUC4 {huc4_id}: {len(resv_in)}")
for _, r in resv_in.iterrows():
    print(f"    {r['dam name']}  ({r['lattitude']:.3f}, {r['longtitude']:.3f})"
          f"  missing={r['missing_ratio']:.3f}")

# ── Step 3: Fetch USGS lake-elevation data for each reservoir ─────────────────
print("\n" + "=" * 55)
print("STEP 3: Fetching USGS reservoir elevation (param 00054)...")
print("=" * 55)

def fetch_elev_for_reservoir(name, lat, lon, radius=0.4):
    """
    Find nearest USGS lake/reservoir site and return daily elevation (ft).
    Tries multiple parameter codes and progressively wider search radii.
    Returns (Series, site_description) or (None, None).
    """
    PARAM_CODES = ["00054", "72020", "62614"]   # reservoir elev, lake stage, lake elev
    SITE_TYPES  = ["LK", "LK,ES", ""]           # lake only, lake+estuary, any

    sites = pd.DataFrame()
    used_param = None

    for r in [radius, radius * 1.5, radius * 2.5]:
        for pcd in PARAM_CODES:
            for stype in SITE_TYPES:
                try:
                    kw = dict(bBox=(lon-r, lat-r, lon+r, lat+r), parameterCd=pcd)
                    if stype:
                        kw["siteType"] = stype
                    result, _ = nwis.get_info(**kw)
                    if not result.empty:
                        sites = result
                        used_param = pcd
                        break
                except Exception:
                    pass
            if not sites.empty:
                break
        if not sites.empty:
            break

    if sites.empty:
        return None, None

    sites["dist"] = (
        (sites["dec_lat_va"].astype(float) - lat) ** 2 +
        (sites["dec_long_va"].astype(float) - lon) ** 2
    ) ** 0.5
    best    = sites.sort_values("dist").iloc[0]
    site_no = best["site_no"]

    try:
        df, _ = nwis.get_dv(
            sites=site_no, parameterCd=used_param,
            start="2000-01-01", end="2020-12-31"
        )
    except Exception:
        return None, None

    if df.empty:
        return None, None

    cols = [c for c in df.columns if used_param in c]
    if not cols:
        cols = [c for c in df.columns if any(p in c for p in PARAM_CODES)]
    if not cols:
        return None, None

    s = pd.to_numeric(df[cols[0]], errors="coerce")
    s.index = pd.to_datetime(s.index)
    return s, f"{site_no} ({best['station_nm']}, param {used_param})"

elev_series = []
for _, row in resv_in.iterrows():
    s, site_desc = fetch_elev_for_reservoir(
        row["dam name"], row["lattitude"], row["longtitude"]
    )
    if s is not None:
        print(f"  ✓ {row['dam name']} → {site_desc}  ({s.notna().sum()} daily values)")
        elev_series.append(s)
    else:
        print(f"  ✗ {row['dam name']} — no USGS lake elevation data found")
    time.sleep(0.3)

# ── Step 4: Compute monthly RDI ───────────────────────────────────────────────
print("\n" + "=" * 55)
print("STEP 4: Computing RDI...")
print("=" * 55)

full_months = pd.date_range("2000-01", "2020-12", freq="MS")

def standardize_per_month(monthly_series, min_hist=2):
    """
    For each month in the series, compute z-score relative to the
    historical values for that calendar month. Returns a Series aligned
    to monthly_series.index.
    """
    rdi_vals = []
    for d in monthly_series.index:
        hist = monthly_series[monthly_series.index.month == d.month].dropna()
        val  = monthly_series.loc[d]
        if pd.isna(val) or len(hist) < min_hist:
            rdi_vals.append(np.nan)
        else:
            mu, sigma = hist.mean(), hist.std()
            rdi_vals.append((val - mu) / sigma if (sigma and sigma > 0) else 0.0)
    return pd.Series(rdi_vals, index=monthly_series.index)

rdi_available = False
rdi = pd.Series(np.nan, index=full_months, name="RDI")

if elev_series:
    # Standardize each reservoir independently, then average
    rdi_parts = []
    for s in elev_series:
        s.index = pd.to_datetime(s.index).tz_localize(None)
        monthly_s = s.resample("MS").mean()
        print(f"    data range: {monthly_s.dropna().index.min().date()} "
              f"→ {monthly_s.dropna().index.max().date()}  "
              f"({monthly_s.notna().sum()} months)")
        zs = standardize_per_month(monthly_s, min_hist=2)
        rdi_parts.append(zs)

    if rdi_parts:
        rdi_raw = pd.concat(rdi_parts, axis=1).mean(axis=1)
        rdi_raw.index = pd.to_datetime(rdi_raw.index).tz_localize(None)
        rdi = rdi_raw.reindex(full_months).rename("RDI")
        rdi_available = rdi.notna().sum() > 12
        print(f"  RDI computed: {rdi.notna().sum()} valid months")

# Fallback: use standardized annual streamflow if RDI still has < 12 months
FLOW_CSV = os.path.join(BASE_DIR, "output", "correlation",
                         "riverflow_resv_HUC4_annual_2000_2020.csv")
if not rdi_available and os.path.exists(FLOW_CSV):
    flow_df = pd.read_csv(FLOW_CSV, dtype={"HUC4": str})
    flow_huc = flow_df[flow_df["HUC4"] == huc4_id].copy()
    if not flow_huc.empty:
        flow_huc["Date"] = pd.to_datetime(flow_huc["Year"].astype(int).astype(str) + "-07-01")
        flow_annual = flow_huc.set_index("Date")["Flow_cfs"].astype(float)
        mu_f, sig_f = flow_annual.mean(), flow_annual.std()
        flow_z = (flow_annual - mu_f) / sig_f if sig_f > 0 else flow_annual * 0
        # Expand annual to monthly (flat within each year)
        rdi_proxy = []
        for d in full_months:
            year_dates = [dt for dt in flow_z.index if dt.year == d.year]
            rdi_proxy.append(flow_z.loc[year_dates[0]] if year_dates else np.nan)
        rdi = pd.Series(rdi_proxy, index=full_months, name="RDI (streamflow proxy)")
        rdi_available = rdi.notna().sum() > 0
        print(f"  Fallback: RDI from annual streamflow Z-score "
              f"({rdi.notna().sum()} months filled)")
    else:
        print(f"  No streamflow data for HUC4 {huc4_id}")
elif not rdi_available:
    print("  No elevation or streamflow data — RDI unavailable")

# ── Step 5: Load SPEI-1 ───────────────────────────────────────────────────────
spei_series = spei[huc4_id].rename("SPEI-1")
print(f"\nSPEI-1 loaded: {spei_series.notna().sum()} valid months")

# ── Step 6: Load WSDI ─────────────────────────────────────────────────────────
wsdi_target = wsdi_df[wsdi_df["HUC4"].str.zfill(4) == huc4_id].copy()
wsdi_target["Date"] = pd.to_datetime(wsdi_target["Year"].astype(int).astype(str) + "-07-01")
wsdi_target = wsdi_target.set_index("Date")["WSDI"].astype(float)
print(f"WSDI loaded: {len(wsdi_target)} annual values")

# ── Step 7: Plot ──────────────────────────────────────────────────────────────
print("\n" + "=" * 55)
print("STEP 7: Plotting...")
print("=" * 55)

fig, axes = plt.subplots(3, 1, figsize=(15, 11), sharex=True,
                          gridspec_kw={"hspace": 0.50})

CLRS = {"spei": "#1565C0", "rdi": "#E65100", "wsdi_pos": "#2E7D32", "wsdi_neg": "#C62828"}
DL   = dict(color="gray", linewidth=0.8, linestyle="--", alpha=0.8)

# ─── Row 1: SPEI-1 ────────────────────────────────────────────────────────────
ax = axes[0]
x  = spei_series.index
y  = spei_series.values
ax.plot(x, y, color=CLRS["spei"], linewidth=1.3, label="SPEI-1")
ax.fill_between(x, y, 0, where=(y >= 0), alpha=0.15, color=CLRS["spei"])
ax.fill_between(x, y, 0, where=(y < 0),  alpha=0.25, color="red")
ax.axhline(0,  color="black", linewidth=0.6)
ax.axhline(-1, **DL, label="Moderate drought (−1.0)")
ax.axhline(1,  **DL)
ax.set_ylabel("SPEI-1", fontsize=10)
ax.set_ylim(-3.5, 3.5)
ax.legend(fontsize=8, loc="lower right")
ax.set_title("SPEI-1  (Standardized Precipitation-Evapotranspiration Index, 1-month)",
             fontsize=9, loc="left", pad=4)
ax.grid(True, linestyle=":", alpha=0.35)

# ─── Row 2: RDI ───────────────────────────────────────────────────────────────
ax = axes[1]
if rdi_available:
    x2 = rdi.index
    y2 = rdi.values
    ax.plot(x2, y2, color=CLRS["rdi"], linewidth=1.3, label="RDI")
    ax.fill_between(x2, y2, 0, where=(np.nan_to_num(y2) >= 0), alpha=0.15, color=CLRS["rdi"])
    ax.fill_between(x2, y2, 0, where=(np.nan_to_num(y2) <  0), alpha=0.25, color="red")
    ax.axhline(0,  color="black", linewidth=0.6)
    ax.axhline(-1, **DL, label="Below-normal storage (−1.0)")
    ax.set_ylim(-3.5, 3.5)
    ax.legend(fontsize=8, loc="lower right")
    resv_names = "; ".join(resv_in["dam name"].tolist())
    ax.set_title(
        f"RDI  (Reservoir Deficit Index — standardized elevation anomaly)\n"
        f"Reservoirs: {resv_names}",
        fontsize=9, loc="left", pad=4
    )
else:
    ax.text(0.5, 0.5,
            "RDI not available: no USGS lake elevation data (param 00054)\n"
            "for reservoirs in this HUC4",
            ha="center", va="center", transform=ax.transAxes,
            fontsize=10, color="gray", style="italic")
    ax.set_title("RDI  (Reservoir Deficit Index — data unavailable)", fontsize=9, loc="left")
ax.set_ylabel("RDI", fontsize=10)
ax.grid(True, linestyle=":", alpha=0.35)

# ─── Row 3: WSDI (annual bars) ────────────────────────────────────────────────
ax = axes[2]
bar_colors = [CLRS["wsdi_pos"] if v >= 0 else CLRS["wsdi_neg"]
              for v in wsdi_target.values]
ax.bar(wsdi_target.index, wsdi_target.values,
       width=270, color=bar_colors, alpha=0.75, align="center")
ax.axhline(0, color="black", linewidth=0.6)
ax.set_ylabel("WSDI", fontsize=10)
ax.set_title(
    f"WSDI  ({WSDI_TYPE.replace('_', ' ')} — annual  (Supply − Demand) / Demand)",
    fontsize=9, loc="left", pad=4
)
legend_handles = [
    mpatches.Patch(color=CLRS["wsdi_pos"], alpha=0.75, label="Supply surplus (WSDI ≥ 0)"),
    mpatches.Patch(color=CLRS["wsdi_neg"], alpha=0.75, label="Supply deficit  (WSDI < 0)"),
]
ax.legend(handles=legend_handles, fontsize=8, loc="lower right")
ax.grid(True, linestyle=":", alpha=0.35)

# ─── Shared x-axis formatting ─────────────────────────────────────────────────
axes[-1].set_xlabel("Year", fontsize=10)
axes[-1].xaxis.set_major_locator(mdates.YearLocator(2))
axes[-1].xaxis.set_major_formatter(mdates.DateFormatter("%Y"))
plt.setp(axes[-1].xaxis.get_majorticklabels(), rotation=30, ha="right")
axes[-1].set_xlim(pd.Timestamp("1999-07-01"), pd.Timestamp("2021-06-01"))

fig.suptitle(
    f"Drought Index Time Series — HUC4 {huc4_id}: {huc4_name}\n"
    f"SPEI-1 (monthly)  |  RDI (monthly)  |  WSDI (annual),  2000–2020",
    fontsize=12, fontweight="bold"
)

out_path = os.path.join(OUT_DIR, f"fig_rdi_wsdi_spei_{huc4_id}.png")
plt.savefig(out_path, dpi=150, bbox_inches="tight")
plt.close()
print(f"\n  Saved → {out_path}")
print("=" * 55)
print("Done.")
