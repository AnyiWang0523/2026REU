"""
Generate Figures 3, 4, 5 for the drought index report.
  Fig 3 — Sample HUC4 monthly time series (rainfall + all 5 indices)
  Fig 4 — PDSI lag analysis (median r vs lag 0-6 months)
  Fig 5 — SPI-1 vs SPEI-1 scatter + difference histogram
"""

import os
import numpy as np
import pandas as pd
import matplotlib.pyplot as plt
import matplotlib.gridspec as gridspec
from scipy import stats

# ── Paths ────────────────────────────────────────────────────────────────────
BASE  = r"C:\Users\anyiw\OneDrive - University of Illinois - Urbana\GIS\2026REU\WSDI"
DATA  = os.path.join(BASE, "output", "drought_correlation")
FIGS  = os.path.join(DATA, "figures")
os.makedirs(FIGS, exist_ok=True)

# ── Load time series ──────────────────────────────────────────────────────────
print("Loading time series data...")
rain  = pd.read_csv(os.path.join(DATA, "huc4_rainfall_mm.csv"),  index_col=0, parse_dates=True)
spi1  = pd.read_csv(os.path.join(DATA, "huc4_spi1.csv"),         index_col=0, parse_dates=True)
spi3  = pd.read_csv(os.path.join(DATA, "huc4_spi3.csv"),         index_col=0, parse_dates=True)
spei1 = pd.read_csv(os.path.join(DATA, "huc4_spei1.csv"),        index_col=0, parse_dates=True)
spei3 = pd.read_csv(os.path.join(DATA, "huc4_spei3.csv"),        index_col=0, parse_dates=True)
pdsi  = pd.read_csv(os.path.join(DATA, "huc4_pdsi.csv"),         index_col=0, parse_dates=True)

# ── Load per-HUC4 correlation summaries ──────────────────────────────────────
corr_spi1  = pd.read_csv(os.path.join(DATA, "corr_spi1.csv"))
corr_spi3  = pd.read_csv(os.path.join(DATA, "corr_spi3.csv"))
corr_spei1 = pd.read_csv(os.path.join(DATA, "corr_spei1.csv"))
corr_spei3 = pd.read_csv(os.path.join(DATA, "corr_spei3.csv"))
corr_pdsi  = pd.read_csv(os.path.join(DATA, "corr_pdsi.csv"))

# Make sure huc4 column is zero-padded string
for df in [corr_spi1, corr_spi3, corr_spei1, corr_spei3, corr_pdsi]:
    df["huc4"] = df["huc4"].astype(str).str.zfill(4)


# ═══════════════════════════════════════════════════════════════════════════════
# FIGURE 3 — Sample HUC4 monthly time series
# ═══════════════════════════════════════════════════════════════════════════════
print("\n-- Figure 3: Sample HUC4 time series --")

# Pick the HUC4 whose SPI-1 r is closest to the overall median (most "typical")
med_r = corr_spi1["r"].median()
idx   = (corr_spi1["r"] - med_r).abs().idxmin()
huc   = corr_spi1.loc[idx, "huc4"]
r_val = corr_spi1.loc[idx, "r"]
print(f"  Selected HUC4: {huc}  (SPI-1 r = {r_val:.3f}, closest to median {med_r:.3f})")

fig, (ax1, ax2) = plt.subplots(2, 1, figsize=(14, 8), sharex=True,
                                gridspec_kw={"height_ratios": [1, 1.6]})

# Upper: bar chart of monthly precipitation
ax1.bar(rain.index, rain[huc], color="steelblue", alpha=0.70, width=25)
ax1.set_ylabel("Precipitation (mm/month)", fontsize=11)
ax1.set_title(f"HUC4 {huc} — Monthly Precipitation and Drought Indices (2000–2020)\n"
              f"[SPI-1 r = {r_val:.2f} vs same-month rainfall]",
              fontsize=12, fontweight="bold")
ax1.yaxis.grid(True, linestyle="--", alpha=0.5)
ax1.set_axisbelow(True)

# Lower: line plot of all five drought indices
styles = [
    (spi1,  "SPI-1",  "green",    1.4, "-"),
    (spei1, "SPEI-1", "darkorange", 1.4, "-"),
    (spi3,  "SPI-3",  "royalblue",  1.1, "--"),
    (spei3, "SPEI-3", "purple",     1.1, "--"),
    (pdsi,  "PDSI",   "crimson",    1.4, "-."),
]
for df, label, color, lw, ls in styles:
    ax2.plot(df.index, df[huc], label=label, color=color,
             linewidth=lw, linestyle=ls, alpha=0.88)

ax2.axhline(0,    color="black", linewidth=0.6, linestyle="-")
ax2.axhline(-0.5, color="gray",  linewidth=0.9, linestyle=":",
            label="Drought threshold (-0.5)")
ax2.set_ylabel("Drought Index Value", fontsize=11)
ax2.set_xlabel("Date", fontsize=11)
ax2.legend(ncol=3, fontsize=9, loc="lower right")
ax2.yaxis.grid(True, linestyle="--", alpha=0.4)
ax2.set_axisbelow(True)

plt.tight_layout()
out3 = os.path.join(FIGS, "timeseries_sample_huc4.png")
plt.savefig(out3, dpi=150, bbox_inches="tight")
plt.close()
print(f"  Saved -> {out3}")


# ═══════════════════════════════════════════════════════════════════════════════
# FIGURE 4 — PDSI lag analysis
# ═══════════════════════════════════════════════════════════════════════════════
print("\n-- Figure 4: PDSI lag analysis --")

# Restrict to 2000-2018 for PDSI (Dai dataset ends Dec 2018)
rain_w  = rain.loc[:"2018-12"]
pdsi_w  = pdsi.loc[:"2018-12"]

# Keep only HUC4s present in both and with <10% NaN in PDSI
common = [h for h in pdsi_w.columns
          if h in rain_w.columns and pdsi_w[h].isna().mean() < 0.10]
print(f"  Valid HUC4s for lag analysis: {len(common)}")

lags      = list(range(0, 7))
med_rs    = []
p25_rs    = []
p75_rs    = []

for lag in lags:
    rs = []
    for h in common:
        p_s = rain_w[h].dropna()
        d_s = pdsi_w[h].dropna()
        # align on common dates, then apply lag: shift rain forward by lag steps
        aligned = pd.DataFrame({"rain": p_s, "pdsi": d_s}).dropna()
        if lag > 0:
            # correlate pdsi at time t with rain at time t-lag
            aligned["rain_lag"] = aligned["rain"].shift(lag)
            aligned = aligned.dropna()
        else:
            aligned["rain_lag"] = aligned["rain"]
        if len(aligned) < 24:
            continue
        r, _ = stats.pearsonr(aligned["rain_lag"], aligned["pdsi"])
        rs.append(r)
    med_rs.append(float(np.median(rs)))
    p25_rs.append(float(np.percentile(rs, 25)))
    p75_rs.append(float(np.percentile(rs, 75)))
    print(f"  Lag {lag}: median r = {med_rs[-1]:.3f}  (n={len(rs)} HUC4s)")

fig, ax = plt.subplots(figsize=(9, 5))
x      = np.array(lags)
colors = ["#2196F3" if i == np.argmax(med_rs) else "#90CAF9" for i in range(len(lags))]
bars   = ax.bar(x, med_rs, color=colors, alpha=0.85, zorder=3)
ax.errorbar(x, med_rs,
            yerr=[np.array(med_rs) - np.array(p25_rs),
                  np.array(p75_rs) - np.array(med_rs)],
            fmt="none", color="navy", capsize=5, linewidth=1.2, zorder=4)

ax.set_xlabel("Precipitation Lag (months prior to PDSI)", fontsize=11)
ax.set_ylabel("Median Pearson r across HUC4s", fontsize=11)
ax.set_title("PDSI vs. Lagged Precipitation — Median Pearson r Across 124 HUC4s (2000–2018)",
             fontsize=11, fontweight="bold")
ax.set_xticks(x)
ax.set_xticklabels(["t (same month)" if k == 0 else f"t-{k}" for k in lags])
ax.axhline(0, color="black", linewidth=0.5)
ax.yaxis.grid(True, linestyle="--", alpha=0.5, zorder=0)
ax.set_axisbelow(True)

# Annotate bar values
for bar, v in zip(bars, med_rs):
    ax.text(bar.get_x() + bar.get_width() / 2, v + 0.005,
            f"{v:.3f}", ha="center", va="bottom", fontsize=8.5)

ax.legend(handles=[
    plt.Rectangle((0, 0), 1, 1, color="#2196F3", alpha=0.85, label="Peak lag"),
    plt.Rectangle((0, 0), 1, 1, color="#90CAF9", alpha=0.85, label="Other lags"),
    plt.Line2D([0], [0], color="navy", linewidth=1.5, label="IQR (25th–75th %ile)"),
], fontsize=9)

plt.tight_layout()
out4 = os.path.join(FIGS, "pdsi_lag_analysis.png")
plt.savefig(out4, dpi=150, bbox_inches="tight")
plt.close()
print(f"  Saved -> {out4}")


# ═══════════════════════════════════════════════════════════════════════════════
# FIGURE 5 — SPI-1 vs SPEI-1 scatter + difference histogram
# ═══════════════════════════════════════════════════════════════════════════════
print("\n-- Figure 5: SPI-1 vs SPEI-1 comparison --")

merged = (corr_spi1[["huc4", "r"]].rename(columns={"r": "r_spi1"})
          .merge(corr_spei1[["huc4", "r"]].rename(columns={"r": "r_spei1"}),
                 on="huc4"))
merged["diff"] = merged["r_spi1"] - merged["r_spei1"]
n_spei_better  = (merged["diff"] < 0).sum()
pct_spei_better = 100 * n_spei_better / len(merged)
med_diff = merged["diff"].median()
print(f"  HUC4s where SPEI-1 > SPI-1: {n_spei_better} ({pct_spei_better:.1f}%)")
print(f"  Median difference (SPI-1 - SPEI-1): {med_diff:.3f}")

fig, (ax1, ax2) = plt.subplots(1, 2, figsize=(12, 5))

# Left: scatter
lim_lo = min(merged[["r_spi1", "r_spei1"]].min()) - 0.05
lim_hi = max(merged[["r_spi1", "r_spei1"]].max()) + 0.05

sc = ax1.scatter(merged["r_spi1"], merged["r_spei1"],
                 c=merged["diff"], cmap="RdYlGn_r",
                 vmin=-0.20, vmax=0.20,
                 s=45, alpha=0.80, edgecolors="white", linewidth=0.3, zorder=3)
ax1.plot([lim_lo, lim_hi], [lim_lo, lim_hi], "k--", linewidth=0.9,
         label="1:1 line (SPI-1 = SPEI-1)")
ax1.set_xlim(lim_lo, lim_hi)
ax1.set_ylim(lim_lo, lim_hi)
ax1.set_xlabel("Pearson r — SPI-1 vs Rainfall", fontsize=11)
ax1.set_ylabel("Pearson r — SPEI-1 vs Rainfall", fontsize=11)
ax1.set_title(f"SPI-1 vs SPEI-1 Correlation per HUC4\n"
              f"(SPEI-1 better in {pct_spei_better:.0f}% of watersheds)",
              fontsize=11, fontweight="bold")
ax1.legend(fontsize=9)
ax1.grid(True, linestyle="--", alpha=0.4)
plt.colorbar(sc, ax=ax1, label="r(SPI-1) - r(SPEI-1)")

# Right: histogram of differences
ax2.hist(merged["diff"], bins=22, color="steelblue",
         alpha=0.75, edgecolor="white", zorder=3)
ax2.axvline(0,        color="black", linestyle="--", linewidth=1.0, label="Zero")
ax2.axvline(med_diff, color="crimson", linestyle="-", linewidth=1.5,
            label=f"Median = {med_diff:.3f}")
ax2.set_xlabel("r(SPI-1) - r(SPEI-1)", fontsize=11)
ax2.set_ylabel("Number of HUC4s", fontsize=11)
ax2.set_title("Distribution of Correlation Differences\n(positive = SPI-1 better)",
              fontsize=11, fontweight="bold")
ax2.legend(fontsize=9)
ax2.grid(True, linestyle="--", alpha=0.4, zorder=0)
ax2.set_axisbelow(True)

plt.tight_layout()
out5 = os.path.join(FIGS, "spi1_spei1_comparison.png")
plt.savefig(out5, dpi=150, bbox_inches="tight")
plt.close()
print(f"  Saved -> {out5}")

print("\nAll three figures generated successfully.")
print(f"  Fig 3 -> timeseries_sample_huc4.png")
print(f"  Fig 4 -> pdsi_lag_analysis.png")
print(f"  Fig 5 -> spi1_spei1_comparison.png")
