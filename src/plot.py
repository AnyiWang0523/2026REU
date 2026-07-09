# ============================================================
# plot.py
# Visualization functions for HUC4 drought-correlation results.
# ============================================================

import matplotlib.pyplot as plt
import pandas as pd
import geopandas as gpd


def plot_correlation_map(huc4_gdf, corr_df, title, ax=None):
    """
    Choropleth map: fill each HUC4 polygon with its Pearson r value.

    Color scale: RdYlGn diverging palette centered at 0.
      Green (r → +1): index tracks rainfall well (positive correlation).
      Red   (r → -1): index moves opposite to rainfall (unexpected).
      Yellow (r ≈ 0): no linear relationship.

    Parameters
    ----------
    huc4_gdf : GeoDataFrame with 'huc4' column and polygon geometry
    corr_df  : pd.DataFrame with index=HUC4 ID and column 'r'
    title    : str, plot title
    ax       : matplotlib Axes (created if None)
    """
    gdf = huc4_gdf.copy()
    gdf = gdf.merge(
        corr_df[["r"]],
        left_on="huc4", right_index=True, how="left"
    )
    if ax is None:
        _, ax = plt.subplots(figsize=(14, 7))

    gdf.plot(
        column="r",
        cmap="RdYlGn",
        vmin=-1, vmax=1,
        legend=True,
        legend_kwds={"label": "Pearson r", "orientation": "horizontal", "shrink": 0.6},
        edgecolor="gray", linewidth=0.3,
        missing_kwds={"color": "lightgrey", "label": "No data"},
        ax=ax,
    )
    ax.set_title(title, fontsize=12)
    ax.set_xlim(-130, -60)
    ax.set_ylim(24, 50)
    ax.axis("off")
    return ax


def plot_correlation_boxplot(summary_df, save_path=None):
    """
    Box plot comparing the distribution of Pearson r values across indices.

    Each box shows the median, IQR, and whiskers (1.5×IQR) of
    per-HUC4 correlation coefficients for one drought index.
    Higher median r → that index tracks rainfall more consistently
    across all US watersheds.

    The horizontal dashed line at r=0 marks 'no correlation'.

    Parameters
    ----------
    summary_df : pd.DataFrame  columns=index names, rows=HUC4 IDs, values=r
    save_path  : str or None, path to save the figure PNG
    """
    fig, ax = plt.subplots(figsize=(10, 5))
    summary_df.boxplot(ax=ax)
    ax.axhline(0, color="gray", linestyle="--", linewidth=0.8)
    ax.set_ylabel("Pearson r  (drought index vs monthly rainfall)")
    ax.set_title("HUC4 Correlation Distribution — Drought Index vs Rainfall\n"
                 "(higher median = index better tracks precipitation variability)")
    plt.tight_layout()
    if save_path:
        plt.savefig(save_path, dpi=150, bbox_inches="tight")
    return fig, ax
