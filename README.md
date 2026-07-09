# Drought Index–Precipitation Correlation Across U.S. HUC4 Watersheds

**2026 REU Project — University of Illinois Urbana-Champaign**

Comparative analysis of three drought indices (SPI, SPEI, PDSI) and their correlation with monthly precipitation across ~222 CONUS HUC4 watersheds (2000–2020).

---

## Research Question

Which drought index most directly tracks same-month precipitation variability at the watershed scale?

**Key Finding:** SPI-1 > SPEI-1 > SPI-3 > SPEI-3 > PDSI (ranked by median Pearson r with same-month rainfall)

---

## Project Structure

```
WSDI/
├── 01_data_prep/
│   └── extract_huc4.py          # Extract HUC4 boundary attributes
├── 02_wsdi/
│   ├── wsdi_ir_gw.py            # Water supply deficit index (irrigation, groundwater)
│   ├── wsdi_ir_sw.py            # Water supply deficit index (irrigation, surface water)
│   ├── wsdi_ps_gw.py            # Water supply deficit index (public supply, groundwater)
│   └── wsdi_ps_sw.py            # Water supply deficit index (public supply, surface water)
├── 03_reservoir_qc/
│   ├── build_reservoir_table.py # Build reservoir coverage table from NID data
│   ├── compute_rdi.py           # Compute Reservoir Demand Index (RDI)
│   └── filter_valid_huc4.py     # Filter HUC4 units with sufficient reservoir coverage
├── 04_drought_analysis/
│   ├── drought_index_analysis.py    # Main correlation analysis (SPI/SPEI/PDSI vs rainfall)
│   └── generate_figures_3_4_5.py    # Generate time series, lag, and SPI vs SPEI figures
├── 05_streamflow/
│   ├── fetch_riverflow_resv.py  # Fetch USGS streamflow data for reservoir watersheds
│   └── correlation_discharge_supply.py  # Correlate discharge with water supply
├── 06_visualization/
│   ├── plot_ir.py               # Irrigation water use plots
│   ├── plot_ps.py               # Public supply water use plots
│   └── plot_rdi_wsdi_spei.py    # Combined RDI/WSDI/SPEI spatial plots
├── src/
│   ├── aggregate.py             # Spatial aggregation (regionmask + OPeNDAP)
│   ├── compute_indices.py       # SPI and SPEI computation from scratch
│   ├── correlation.py           # Pearson r and lag correlation utilities
│   └── plot.py                  # Shared plotting helpers
├── drought_index_notebook.ipynb     # Main analysis notebook (full pipeline)
└── reservoir_storage_comparison.ipynb
```

---

## Data Sources

| Dataset | Source | Resolution | Period |
|---|---|---|---|
| Monthly Precipitation (CMAP Enhanced) | NOAA PSL OPeNDAP | 2.5°, monthly | 2000–2020 |
| Monthly Temperature (GHCN+CAMS) | NOAA PSL OPeNDAP | 0.5°, monthly | 2000–2020 |
| PDSI (Dai scPDSI) | NOAA PSL OPeNDAP | 2.5°, monthly | 2000–2018 |
| HUC4 Boundaries | USGS NHD WBD | — | — |
| National Inventory of Dams | USACE NID | — | — |

All climate data accessed remotely via OPeNDAP — no large local downloads required.

> **Note:** Raw data files (`data/`) and generated outputs (`output/`) are not tracked in this repository due to size. See the data sources above to reproduce them.

---

## Methods

### Drought Indices

- **SPI** (McKee et al. 1993): Gamma distribution fit to rolling precipitation accumulation (1-month, 3-month), transformed to standard normal via probit.
- **SPEI** (Vicente-Serrano et al. 2010): Same as SPI but applied to the water balance D = P − PET. PET computed via Thornthwaite (1948) using monthly temperature and latitude.
- **PDSI** (Palmer 1965, Dai self-calibrated): Pre-computed product from NOAA PSL. Has ~10-month time constant due to recursive smoothing.

### Spatial Aggregation

HUC4 polygon masks built with `regionmask`. Each grid cell assigned to the HUC4 it falls within; per-basin time series computed as arithmetic mean.

### Correlation Analysis

Per-HUC4 Pearson r between each drought index and same-month precipitation. Two-sided significance test (α = 0.05). PDSI lag analysis: precipitation shifted 0–6 months back to characterize PDSI's memory effect.

---

## Key Results

| Index | Median Pearson r | % Significant HUC4s |
|---|---|---|
| SPI-1 | ~0.76 | ~97% |
| SPEI-1 | ~0.69 | ~95% |
| SPI-3 | ~0.55 | ~90% |
| SPEI-3 | ~0.50 | ~88% |
| PDSI | ~0.34 | ~72% |

- **PDSI** peak correlation occurs at ~2–3 month lag, confirming its low-pass filter behavior.
- **SPI-1** is theoretically a monotonic transform of same-month precipitation, so its near-1.0 Spearman ρ is expected.
- **SPEI** is preferred for operational monitoring because it accounts for evaporative demand — a no-rain winter month with low PET is not a drought.

---

## Setup

```bash
python -m venv .venv
.venv\Scripts\activate        # Windows
pip install numpy pandas xarray scipy regionmask matplotlib geopandas pydap
```

No additional data downloads needed for the OPeNDAP-based analysis. Local NID and HUC4 shapefiles are required for `03_reservoir_qc/` and `01_data_prep/`.

---

## References

- McKee, T.B., Doesken, N.J., & Kleist, J. (1993). *The relationship of drought frequency and duration to time scales.* Proc. 8th Conference on Applied Climatology, AMS.
- Palmer, W.C. (1965). *Meteorological Drought.* U.S. Weather Bureau Research Paper No. 45.
- Vicente-Serrano, S.M., Beguería, S., & López-Moreno, J.I. (2010). A multiscalar drought index sensitive to global warming: The SPEI. *Journal of Climate, 23*(7), 1696–1718.
- Dai, A. (2011). Characteristics and trends in various forms of the PDSI during 1900–2008. *Journal of Geophysical Research: Atmospheres, 116*, D12115.
- Thornthwaite, C.W. (1948). An approach toward a rational classification of climate. *Geographical Review, 38*(1), 55–94.
