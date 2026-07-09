# ============================================================
# correlation.py
# Per-HUC4 time-series correlation between drought indices
# and monthly rainfall.
# ============================================================

import numpy as np
import pandas as pd
from scipy import stats


def compute_correlation(df_index, df_rain, method="pearson"):
    """
    Compute per-HUC4 correlation between a drought index and rainfall.

    Two methods are offered:

    Pearson r  (linear correlation):
    ---------------------------------
      r = Σ[(xᵢ - x̄)(yᵢ - ȳ)]
          ─────────────────────────────────────
          sqrt[Σ(xᵢ - x̄)²] · sqrt[Σ(yᵢ - ȳ)²]

      Ranges from -1 (perfect inverse) to +1 (perfect direct).
      Measures how well the relationship fits a straight line.
      Assumes approximate normality; sensitive to outliers.

    Spearman ρ  (rank correlation):
    ---------------------------------
      ρ = 1 - 6·Σdᵢ²
              ───────────
              n(n² - 1)
      where dᵢ = rank(xᵢ) - rank(yᵢ).

      Measures monotonic (not necessarily linear) agreement.
      Equivalent to Pearson r computed on the ranks of x and y.
      More robust to outliers and skewed distributions.

    Expected sign:
      SPI and SPEI are standardized so that wet months → positive
      values and dry months → negative.  We therefore expect r > 0
      (more rain → higher index).
      PDSI may lag rainfall by 1–3 months due to its soil-moisture
      'memory' (recursive 0.897 weighting), so its r may be slightly
      lower than SPI-1/SPEI-1 when correlated with same-month rainfall.

    Parameters
    ----------
    df_index : pd.DataFrame  rows=time, cols=HUC4 IDs  (drought index)
    df_rain  : pd.DataFrame  rows=time, cols=HUC4 IDs  (mm/month)
    method   : 'pearson' or 'spearman'

    Returns
    -------
    pd.DataFrame  index=HUC4 ID, columns=['r', 'p', 'n']
      r  : correlation coefficient
      p  : two-sided p-value (H₀: r = 0)
      n  : number of valid paired observations used
    """
    common_idx = df_index.index.intersection(df_rain.index)
    common_huc = df_index.columns.intersection(df_rain.columns)
    df_index = df_index.loc[common_idx, common_huc]
    df_rain  = df_rain.loc[common_idx, common_huc]

    corr_fn = stats.pearsonr if method == "pearson" else stats.spearmanr
    records = []

    for huc in common_huc:
        x = df_index[huc]
        y = df_rain[huc]

        # Drop any time steps where either series is NaN
        # (first `scale-1` months of SPI/SPEI are NaN by construction)
        valid = x.notna() & y.notna()
        x_v, y_v = x[valid].values, y[valid].values

        if len(x_v) < 10:
            records.append({"huc4": huc, "r": np.nan, "p": np.nan, "n": len(x_v)})
            continue

        r, p = corr_fn(x_v, y_v)
        records.append({"huc4": huc, "r": float(r), "p": float(p), "n": int(len(x_v))})

    return pd.DataFrame(records).set_index("huc4")
