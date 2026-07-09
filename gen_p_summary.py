import pandas as pd
import numpy as np
import os
import geopandas as gpd
import matplotlib.pyplot as plt
import matplotlib.patches as mpatches

BASE_DIR = r"C:\Users\anyiw\OneDrive - University of Illinois - Urbana\GIS\2026REU\WSDI"
OUT_DIR  = os.path.join(BASE_DIR, "output", "drought_correlation")
FIG_DIR  = os.path.join(OUT_DIR, "figures")
HUC4_SHP = os.path.join(BASE_DIR, "data", "HUC4", "WBDHU4.shp")
os.makedirs(FIG_DIR, exist_ok=True)

indices = {
    "SPI-1":  "spi1",
    "SPI-3":  "spi3",
    "SPEI-1": "spei1",
    "SPEI-3": "spei3",
    "PDSI":   "pdsi",
}

corr_results = {}
for name, safe in indices.items():
    df = pd.read_csv(os.path.join(OUT_DIR, "corr_" + safe + ".csv"), index_col=0)
    corr_results[name] = df

# ── p-value summary CSV ───────────────────────────────────────────────────────
p_summary = pd.DataFrame({name: corr["p"] for name, corr in corr_results.items()})
p_summary.index.name = "huc4"
p_summary.to_csv(os.path.join(OUT_DIR, "correlation_p_summary.csv"))
print("Saved: correlation_p_summary.csv  shape =", p_summary.shape)

# ── Stats table (Table 2) ─────────────────────────────────────────────────────
rows = []
for name, corr in corr_results.items():
    r = corr["r"].dropna()
    p = corr["p"].dropna()
    n_total = len(r)
    n_sig   = int((p < 0.05).sum())
    rows.append({
        "Index":           name,
        "Median_r":        round(r.median(), 3),
        "Mean_r":          round(r.mean(),   3),
        "Std_r":           round(r.std(),    3),
        "N_HUC4":          n_total,
        "N_Significant":   n_sig,
        "Pct_Significant": round(100.0 * n_sig / n_total, 1),
    })

stats_table = pd.DataFrame(rows)
stats_table.to_csv(os.path.join(OUT_DIR, "stats_table.csv"), index=False)
print("\nStats table (Table 2):")
print(stats_table.to_string(index=False))

# ── Significance map ──────────────────────────────────────────────────────────
huc4 = gpd.read_file(HUC4_SHP).to_crs("EPSG:4326")
huc4["huc4"] = huc4["huc4"].astype(str).str.zfill(4)
huc4 = huc4[huc4["huc4"].str[:2].astype(int) <= 18].copy()

fig, axes = plt.subplots(2, 3, figsize=(22, 12))
axes = axes.flatten()

for i, (name, corr) in enumerate(corr_results.items()):
    ax = axes[i]
    p_df = corr[["p"]].reset_index()
    p_df["huc4"] = p_df["huc4"].astype(str).str.zfill(4)
    merged = huc4.merge(p_df, on="huc4", how="left")

    no_data  = merged[merged["p"].isna()]
    sig      = merged[merged["p"] < 0.05]
    not_sig  = merged[(merged["p"] >= 0.05) & merged["p"].notna()]

    if not no_data.empty:  no_data.plot( ax=ax, color="#BDBDBD", linewidth=0.3)
    if not sig.empty:     sig.plot(     ax=ax, color="#2196F3", linewidth=0.3)
    if not not_sig.empty: not_sig.plot( ax=ax, color="#F44336", linewidth=0.3)

    n_sig   = len(sig)
    n_analyzed = len(sig) + len(not_sig)
    ax.set_title(
        f"{name}  —  p<0.05: {n_sig}/{n_analyzed} analyzed HUC4s"
        f"  ({100*n_sig/n_analyzed:.0f}%)",
        fontsize=9
    )
    ax.axis("off")
    ax.legend(handles=[
        mpatches.Patch(color="#2196F3", label="p < 0.05 (significant)"),
        mpatches.Patch(color="#F44336", label="p >= 0.05 (not significant)"),
        mpatches.Patch(color="#BDBDBD", label="Not analyzed"),
    ], loc="lower left", fontsize=7)

axes[-1].set_visible(False)
plt.suptitle(
    "Pearson r Significance Map — Drought Indices vs Same-Month Rainfall (alpha=0.05)",
    fontsize=12, fontweight="bold"
)
plt.tight_layout()
out_sig = os.path.join(FIG_DIR, "map_significance_all.png")
plt.savefig(out_sig, dpi=150, bbox_inches="tight")
plt.close()
print("\nSaved: map_significance_all.png")
print("Done.")
