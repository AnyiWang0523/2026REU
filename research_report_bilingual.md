# Comparative Assessment of Drought Index–Precipitation Correlation Across U.S. HUC4 Watersheds
# 美国HUC4流域干旱指数与降水相关性比较评估

**Running Head / 页眉缩写:** DROUGHT INDEX–PRECIPITATION CORRELATION IN HUC4 WATERSHEDS

**Author / 作者:** [Your Name]

**Institutional Affiliation / 单位:**
Department of Geography & Geographic Information Science,
University of Illinois Urbana-Champaign
伊利诺伊大学厄巴纳-香槟分校 地理与地理信息科学系

**Program / 项目:** 2026 Research Experience for Undergraduates (REU) — University of Illinois Urbana-Champaign

**Author Note / 作者说明:**
This research was conducted as part of the 2026 REU program at the University of Illinois Urbana-Champaign. The author thanks the program advisors and the open-data repositories maintained by NOAA PSL and USGS for making this work possible. Correspondence: anyiw887@gmail.com.

本研究在伊利诺伊大学厄巴纳-香槟分校2026年本科生科研体验（REU）项目框架下完成。作者感谢项目导师以及NOAA PSL和USGS维护的开放数据库对本研究的支持。通讯邮件：anyiw887@gmail.com。

---

## Abstract / 摘要

Selecting an appropriate drought index is critical for operational water resource management and sub-seasonal to seasonal (S2S) hydrological forecasting, yet systematic comparisons across climatologically diverse watersheds remain limited. This study evaluates how well five commonly used drought indices — the Standardized Precipitation Index at 1- and 3-month scales (SPI-1, SPI-3), the Standardized Precipitation-Evapotranspiration Index at the same scales (SPEI-1, SPEI-3), and the Palmer Drought Severity Index (PDSI) — track concurrent monthly precipitation across 124 CONUS HUC4 watersheds with valid multi-index data coverage over the period 2000–2020. All drought indices were computed from remotely accessed NOAA PSL gridded climate products and spatially aggregated to HUC4 boundaries using polygon masking. Pearson correlation coefficients were calculated per watershed between each index and same-month precipitation, and a lag analysis (0–6 months) was conducted for PDSI. Results show a clear ranking: SPI-1 (median r ≈ 0.76) > SPEI-1 (≈ 0.69) > SPI-3 (≈ 0.48) > SPEI-3 (≈ 0.45) > PDSI (≈ 0.34). PDSI's same-month correlation is the highest among its own lags (r = 0.34 at lag 0), but declines monotonically to r ≈ 0.21 at a 6-month lag, reflecting the progressive dilution of the precipitation signal by the recursive memory structure. These findings provide evidence-based guidance for index selection: SPI-1 is the most direct precipitation proxy; SPEI-1 is preferred for operational drought monitoring; and PDSI best represents cumulative soil moisture deficits.

*Keywords: drought index, SPI, SPEI, PDSI, HUC4, precipitation correlation, watershed hydrology*

---

选择合适的干旱指数对于水资源运营管理和次季节至季节（S2S）水文预报至关重要，但针对气候多样流域的系统性比较研究仍然有限。本研究评估了五种常用干旱指数——1个月和3个月尺度的标准化降水指数（SPI-1、SPI-3）、同尺度的标准化降水蒸散指数（SPEI-1、SPEI-3）以及帕尔默干旱烈度指数（PDSI）——在2000–2020年间与124个具备有效多指数数据覆盖的CONUS HUC4流域同月降水的同步追踪能力。所有干旱指数均由远程访问的NOAA PSL格网气候产品计算，并通过多边形掩膜空间聚合至HUC4流域边界。对每个流域计算各指数与同月降水的Pearson相关系数，并对PDSI进行0–6个月滞后分析。结果显示相关性排序为：SPI-1（中位数r ≈ 0.76）> SPEI-1（≈ 0.69）> SPI-3（≈ 0.48）> SPEI-3（≈ 0.45）> PDSI（≈ 0.34）。PDSI的同月相关在自身各滞后中最高（lag 0时r = 0.34），但随滞后时步增加单调递减至约r ≈ 0.21（6个月滞后），反映了递归记忆结构对降水信号的逐步稀释。本研究为指数选择提供了基于证据的依据：SPI-1是最直接的降水代理指标；SPEI-1更适合业务化干旱监测；PDSI最能表征累积土壤水分亏缺。

*关键词：干旱指数、SPI、SPEI、PDSI、HUC4、降水相关性、流域水文学*

---

## 1. Introduction / 引言

### 1.1 Background and Motivation / 研究背景与动机

Drought is one of the most economically and ecologically damaging natural hazards in the United States, imposing cumulative costs exceeding $250 billion since 1980 (NOAA NCEI, 2023). Its impacts cascade from sustained precipitation deficits and elevated evaporative demand into reservoir storage depletion, irrigation shortfalls, reduced hydropower generation, and degraded aquatic ecosystems. Effective drought management therefore depends critically on timely and accurate drought monitoring — the ability to characterize the current state of water availability at the watershed scale and to anticipate how that state will evolve.

干旱是美国最具经济和生态破坏力的自然灾害之一，自1980年以来造成的累计损失超过2500亿美元（NOAA NCEI, 2023）。持续的降水亏缺和升高的蒸发需求可引发水库蓄水减少、灌溉短缺、水力发电减少以及水生生态系统退化等一系列连锁影响。因此，有效的干旱管理在很大程度上依赖于及时准确的干旱监测——即在流域尺度上表征当前水资源可用状态并预判其演变趋势的能力。

Recent decades have seen substantial advances in both reservoir operation modeling and sub-seasonal to seasonal (S2S) meteorological forecasting. Ensemble streamflow prediction systems, data-driven inflow models, and machine learning approaches have improved the skill of seasonal outlooks for reservoir operations (Yuan et al., 2015). Similarly, dynamical and statistical S2S forecast products now provide skillful precipitation and temperature outlooks at 2–6-week lead times over much of CONUS (Robertson et al., 2015). The practical value of these advances depends, however, on whether the drought characterization inputs used to initialize or validate models accurately reflect the physical conditions that constrain water availability.

近年来，水库运营模型和次季节至季节（S2S）气象预报均取得了显著进展。集合径流预报系统、数据驱动入流模型和机器学习方法提高了水库运营季节预报的技巧（Yuan et al., 2015）。与此同时，动力学和统计S2S预报产品现已能在美国本土大部分地区提供2–6周预见期内技巧较高的降水和气温预报（Robertson et al., 2015）。然而，上述进展的实际应用价值，取决于初始化或验证模型所采用的干旱表征输入是否准确反映了制约水资源可用性的物理条件。

A prerequisite for linking S2S forecasts to water resource management is identifying which drought index most faithfully represents the precipitation variability that drives inflow anomalies. Dozens of indices exist in the literature, each with different temporal integration windows, physical assumptions, and data requirements. The three most widely used are: (1) the Standardized Precipitation Index (SPI; McKee et al., 1993); (2) the Standardized Precipitation-Evapotranspiration Index (SPEI; Vicente-Serrano et al., 2010); and (3) the Palmer Drought Severity Index (PDSI; Palmer, 1965).

将S2S预报与水资源管理相结合的前提，是确定哪种干旱指数最能真实反映驱动入流异常的降水变率。文献中存在数十种干旱指数，各自具有不同的时间积分窗口、物理假设和数据需求。应用最广泛的三种为：(1) 标准化降水指数（SPI；McKee et al., 1993）；(2) 标准化降水蒸散指数（SPEI；Vicente-Serrano et al., 2010）；以及(3) 帕尔默干旱烈度指数（PDSI；Palmer, 1965）。

### 1.2 Research Objectives / 研究目标

This study addresses the following core question: *Which drought index most immediately tracks concurrent monthly precipitation variability at the HUC4 watershed scale across CONUS?*

本研究的核心问题为：*在CONUS HUC4流域尺度上，哪种干旱指数与同月降水变率的同步相关性最高？*

Specific objectives are:
具体目标包括：

1. To compute five drought index time series (SPI-1, SPI-3, SPEI-1, SPEI-3, PDSI) for 124 CONUS HUC4 watersheds with valid coverage across all five indices over 2000–2020 from publicly accessible gridded climate datasets.
   基于公开格网气候数据集，计算2000–2020年间124个CONUS HUC4流域（具备所有五种指数有效覆盖的流域）的五种干旱指数时间序列（SPI-1、SPI-3、SPEI-1、SPEI-3、PDSI）。

2. To quantify the Pearson correlation between each index and same-month precipitation per HUC4 and summarize the spatial distribution of these correlations.
   逐流域量化各指数与同月降水的Pearson相关系数，并汇总相关性的空间分布特征。

3. To characterize the temporal lag structure of the PDSI–precipitation relationship.
   表征PDSI与降水关系的时间滞后结构。

4. To compare SPI and SPEI at matched accumulation scales and identify regions where they diverge.
   在相同积累尺度下比较SPI和SPEI，识别二者相关性差异显著的区域。

### 1.3 Significance / 研究意义

The results provide an evidence-based framework for selecting drought index inputs for reservoir inflow forecasting models, and contribute a drought characterization framework for hydrological model calibration across climatologically diverse U.S. basins.

本研究结果为水库入流预报模型的干旱指数输入选择提供了基于证据的框架，并为跨气候多样性美国流域的水文模型率定提供了干旱表征框架。

---

## 2. Data and Study Area / 数据与研究区域

### 2.1 Study Area / 研究区域

The study area encompasses CONUS HUC4 watersheds as defined by the USGS Watershed Boundary Dataset (WBD). CONUS contains approximately 222 total HUC4 polygons; after requiring valid grid-cell coverage across all five drought indices (constrained primarily by the 2.5° PDSI grid), 124 watersheds were retained in this analysis. These span a wide range of climatic conditions — from the humid Eastern Seaboard (HUC2 regions 01–03) to the semi-arid Great Plains (regions 10–11) and the arid Southwest (regions 14–16). HUC4 boundaries were obtained from the USGS National Hydrography Dataset as a shapefile (EPSG:4326).

研究区域涵盖USGS流域边界数据集（WBD）所定义的CONUS HUC4流域。CONUS共有约222个HUC4多边形；在要求五种干旱指数均具有有效格网覆盖（主要受2.5° PDSI格网限制）后，本分析共保留124个流域。这些流域涵盖从湿润的东部海岸（HUC2区域01–03）到半干旱的大平原（区域10–11）以及干旱的西南地区（区域14–16）等广泛气候条件。HUC4边界由USGS国家水文数据集获取（坐标系EPSG:4326）。

### 2.2 Climate Data / 气候数据

All climate data were accessed remotely via OPeNDAP from NOAA Physical Sciences Laboratory (PSL) servers. Table 1 summarizes the three datasets used.

所有气候数据均通过OPeNDAP协议从NOAA物理科学实验室（PSL）服务器远程获取。表1汇总了三个数据集的基本信息。

**Table 1 / 表1.** Climate datasets used in this study. / 本研究使用的气候数据集。

| Dataset / 数据集 | Source / 来源 | Resolution / 分辨率 | Period / 时段 | Variable / 变量 |
|---|---|---|---|---|
| CMAP Enhanced Precipitation / 增强降水数据 | NOAA PSL OPeNDAP | 2.5°, monthly / 月 | 2000–2020 | Precipitation mm/day / 降水量 |
| GHCN+CAMS Temperature / 气温数据 | NOAA PSL OPeNDAP | 0.5°, monthly / 月 | 2000–2020 | Mean air temperature °C / 月均气温 |
| Dai scPDSI | NOAA PSL OPeNDAP | 2.5°, monthly / 月 | 2000–2018 | Self-calibrated PDSI / 自校准PDSI |

**CMAP Enhanced Precipitation** (Xie & Arkin, 1997): Monthly precipitation estimates merging rain gauge observations, satellite retrievals, and reanalysis model output. Original units (mm/day) were converted to mm/month by multiplying by the number of days in each calendar month.

**CMAP增强降水数据**（Xie & Arkin, 1997）：融合雨量站观测、卫星反演和再分析模式输出的月降水估算产品。原始单位为mm/day，乘以各历月天数转换为mm/month。

**GHCN+CAMS Temperature**: Gridded monthly mean temperature combining Global Historical Climatology Network records with the Climate Anomaly Monitoring System. Used exclusively for computing Thornthwaite PET in the SPEI pipeline.

**GHCN+CAMS气温数据**：融合全球历史气候网络站点记录与气候异常监测系统的格网月均气温产品，仅用于SPEI计算流程中的Thornthwaite PET估算。

**Dai scPDSI** (Dai, 2011): Self-calibrating PDSI, used directly without recomputation. Coverage ends December 2018.

**Dai自校准PDSI**（Dai, 2011）：直接使用，不重新计算，数据截至2018年12月。

### 2.3 Notes on Analysis Period / 分析时段说明

The primary analysis window is January 2000 through December 2020 for SPI and SPEI, aligned with existing WSDI water-use records. For PDSI, the window is truncated at December 2018 due to dataset availability, yielding an 18-year analysis period for PDSI correlations.

SPI和SPEI的主要分析窗口为2000年1月至2020年12月，与现有WSDI用水记录对齐。受数据可用性限制，PDSI分析窗口截至2018年12月，对应18年的PDSI相关性分析时段。

---

## 3. Methods / 方法

### 3.1 Data Retrieval and Spatial Aggregation / 数据获取与空间聚合

All netCDF datasets were retrieved via the `pydap` engine within `xarray`, which avoids compatibility issues between the Windows platform and the default `libnetcdf`/`libcurl` backend. Following retrieval, two coordinate normalization steps were applied universally:

所有netCDF数据集均通过`xarray`的`pydap`引擎获取，以避免Windows平台与默认`libnetcdf`/`libcurl`后端之间的兼容性问题。获取后统一执行两步坐标标准化：

1. **Longitude conversion / 经度转换**: Coordinates on the 0–360° convention were converted to −180–180° via `lon_conv = ((lon + 180) % 360) - 180` and re-sorted ascending. / 将0–360°约定的经度坐标通过`lon_conv = ((lon+180)%360)-180`转换为−180–180°并重新升序排列。

2. **Latitude sorting / 纬度排序**: Latitude arrays were normalized to ascending (south-to-north) order for consistent slicing. / 将纬度数组统一规范为升序（南→北）排列，确保切片一致性。

Spatial aggregation of gridded fields to HUC4 polygons was performed using the `regionmask` library. A boolean mask array of shape (n_HUC4, n_lat, n_lon) was constructed, and the per-HUC4 time series was computed as the arithmetic mean of all assigned grid cells.

格网数据向HUC4多边形的空间聚合采用`regionmask`库实现。构建形状为(n_HUC4, n_lat, n_lon)的布尔掩膜数组，并将分配至各流域的格网单元进行算术均值聚合，得到逐流域时间序列。

### 3.2 Drought Index Computation / 干旱指数计算

**Standardized Precipitation Index (SPI) / 标准化降水指数**

SPI was computed following McKee et al. (1993) at 1-month (SPI-1) and 3-month (SPI-3) accumulation scales:

SPI遵循McKee et al.（1993）方法，按1个月（SPI-1）和3个月（SPI-3）积累尺度计算：

1. *Accumulation / 积累*: Monthly precipitation accumulated over a rolling 1- or 3-month window. / 对月降水按滑动1个月或3个月窗口进行累加。
2. *Distribution fitting / 分布拟合*: Per calendar month, a two-parameter Gamma distribution was fitted to the calibration period (2000–2020) accumulated values, with a mixed distribution to accommodate zero precipitation. / 逐历月对校准期（2000–2020）积累降水拟合双参数Gamma分布，引入零降水混合分布处理无降水月份。
3. *Probit transform / Probit变换*: SPI = Φ⁻¹(H(x)), where H(x) = p₀ + (1 − p₀) · F_Gamma(x). / 通过累积概率的正态分位数变换得到SPI值。

**Standardized Precipitation-Evapotranspiration Index (SPEI) / 标准化降水蒸散指数**

SPEI was computed following Vicente-Serrano et al. (2010) at 1- and 3-month scales:

SPEI遵循Vicente-Serrano et al.（2010）方法，按1个月和3个月尺度计算：

1. *PET computation / PET计算*: Monthly potential evapotranspiration computed by the Thornthwaite (1948) method using monthly mean temperature and latitude. / 采用Thornthwaite（1948）方法，基于月均气温和纬度计算月潜在蒸散量（PET）。
2. *Water balance / 水量平衡*: D = P − PET (mm/month). / 计算月水量平衡：D = P − PET（mm/月）。
3. *Accumulation / 积累*: D accumulated over a rolling 1- or 3-month window. / 对D按滑动1个月或3个月窗口进行累加。
4. *Distribution fitting and transform / 分布拟合与变换*: Accumulated D fitted to a Log-Logistic distribution (approximated by Pearson Type III), then probit-transformed to standard normal quantiles. / 对累积D拟合Log-Logistic分布（由Pearson III近似），再做Probit变换得到SPEI值。

**Palmer Drought Severity Index (PDSI) / 帕尔默干旱烈度指数**

The Dai scPDSI product (Dai, 2011) was used directly without recomputation. The original PDSI update equation is:

直接使用Dai自校准PDSI产品（Dai, 2011），不重新计算。原始PDSI递归更新方程为：

> PDSI_t = 0.897 · PDSI_{t−1} + Z_t / 3

The autoregressive coefficient 0.897 produces a theoretical time constant of approximately 9.5 months, meaning only ~10% of the current PDSI reflects same-month precipitation conditions.

自回归系数0.897产生约9.5个月的理论时间常数，意味着当月PDSI仅约10%反映同月降水状况。

### 3.3 Correlation Analysis / 相关性分析

For each HUC4 and each drought index, Pearson r was computed between the index time series and same-month precipitation over the analysis window (2000–2020 for SPI/SPEI; 2000–2018 for PDSI). Two-sided p-values were computed (α = 0.05). A PDSI lag analysis was performed by shifting precipitation backward by 0–6 months and computing median r across all HUC4s at each lag.

 
---

## 4. Results / 结果

### 4.1 Summary Statistics of Index–Precipitation Correlations / 各指数与同月降水相关性汇总

Table 2 presents per-HUC4 Pearson r statistics for all five drought indices.

表2列出了全部五种干旱指数的逐HUC4 Pearson r统计结果。

**Table 2 / 表2.** Summary of Pearson r between drought indices and same-month precipitation (n = 124 CONUS HUC4, 2000–2020). / 干旱指数与同月降水Pearson r汇总（n = 124个CONUS HUC4，2000–2020）。

| Index / 指数 | Median r / 中位数r | Mean r / 均值r | Std. Dev. / 标准差 | % Significant / 显著比例 (p < 0.05) |
|---|---|---|---|---|
| SPI-1 | 0.76 | 0.74 | 0.13 | 100% |
| SPEI-1 | 0.69 | 0.69 | 0.13 | 100% |
| SPI-3 | 0.48 | 0.48 | 0.08 | 100% |
| SPEI-3 | 0.45 | 0.44 | 0.08 | 100% |
| PDSI | 0.34 | 0.32 | 0.13 | 77% |

The ranking from highest to lowest same-month correlation is: SPI-1 > SPEI-1 > SPI-3 > SPEI-3 > PDSI. SPI-1 (median r = 0.76) and SPEI-1 (0.69) achieve 100% significant watersheds. The 3-month indices are notably lower (SPI-3: 0.48; SPEI-3: 0.45) but also 100% significant. PDSI has the lowest median r (0.34) and the highest proportion of non-significant watersheds (~23%).

同月相关性由高到低的排序为：SPI-1 > SPEI-1 > SPI-3 > SPEI-3 > PDSI。SPI-1（中位r = 0.76）和SPEI-1（0.69）在全部流域均达到显著性。3个月指数明显偏低（SPI-3：0.48；SPEI-3：0.45），但同样在全部流域显著。PDSI的中位r最低（0.34），不显著流域的比例最高（约23%）。

### 4.2 Spatial Distribution of Correlations / 相关系数空间分布（Figure 1 / 图1）

Choropleth maps (Figure 1) reveal that for SPI-1 and SPEI-1, correlations are highest (r > 0.80) in the semi-arid Southwest (HUC2 regions 14–16) and the Southern Plains (region 11), where precipitation is highly episodic and individual monthly totals dominate the drought signal. Correlations are somewhat lower (r = 0.60–0.75) in the humid Southeast and Northeast, where precipitation is more uniformly distributed throughout the year.

分级色彩图（图1）显示，SPI-1和SPEI-1在半干旱的西南地区（HUC2区域14–16）和南部平原（区域11）相关性最高（r > 0.80），这些地区降水具有高度间歇性，单月降水量主导干旱信号。在降水全年分布较均匀的湿润东南部和东北部，相关性略低（r = 0.60–0.75）。

For PDSI, the spatial pattern is more heterogeneous, with low median r across most of CONUS and near-zero or negative r in portions of the Pacific Northwest and high-elevation Rocky Mountain watersheds, where snowmelt dynamics decouple precipitation timing from soil moisture signals.

PDSI的空间格局更为异质，在美国本土大部分地区中位r较低，在太平洋西北部和落基山脉高海拔流域部分区域，r接近于零甚至为负值，这与积雪融水动态将降水时序与土壤水分信号解耦有关。

### 4.3 Distribution of Pearson r Across HUC4s / 各HUC4相关系数分布（Figure 2 / 图2）

Box plots (Figure 2) show that the 3-month indices (SPI-3, SPEI-3) have the most consistent performance across watersheds (std ≈ 0.08), narrower than both the 1-month indices and PDSI (std ≈ 0.13 each). SPI-1 and SPEI-1 achieve higher medians but with wider spread, reflecting greater sensitivity to climate regime differences. PDSI has both the lowest median and notable outliers in the negative range, corresponding to snow-dominated Pacific Northwest basins where PDSI can be weakly or inversely correlated with concurrent precipitation.

箱型图（图2）表明，3个月指数（SPI-3、SPEI-3）在各流域间的表现最为稳定（std ≈ 0.08），离散程度小于1个月指数和PDSI（各约0.13）。SPI-1和SPEI-1中位数更高，但离散程度也更大，反映其对气候区差异的较高敏感性。PDSI中位数最低，且存在显著的负值异常点，对应太平洋西北地区以积雪融水主导的流域——在这些流域中，PDSI与同月降水的相关性可能趋近于零甚至为负。

### 4.4 Sample Time Series / 示例流域时间序列（Figure 3 / 图3）

A representative Great Lakes HUC4 watershed (HUC4 0424, SPI-1 r = 0.76, closest to the median) illustrates the temporal correspondence over 2000–2020 (Figure 3). SPI-1 closely mirrors the precipitation trace with minimal lag. SPEI-1 follows a similar pattern but shows attenuated peaks in summer months when high PET partially offsets precipitation. PDSI is notably smoother and slower-moving: extreme drought events (e.g., 2002–2003, where PDSI reached approximately −8) persist for many months beyond the precipitation deficit that triggered them, and the recovery is gradual rather than coincident with the return of precipitation.

大湖区代表性HUC4流域（HUC4 0424，SPI-1 r = 0.76，最接近中位数）的2000–2020年时间序列（图3）表明：SPI-1与降水序列几乎同步吻合；SPEI-1表现相似，但夏季月份峰值因高PET抵消降水而有所衰减；PDSI明显更平滑、变化更缓慢——极端干旱事件（如2002–2003年PDSI约达−8）在引发降水亏缺后持续数月，且恢复过程是渐进的，而非随降水恢复而同步消退。

### 4.5 PDSI Lag Analysis / PDSI滞后分析（Figure 4 / 图4）

The median Pearson r between PDSI and lagged precipitation (Figure 4) is highest at lag 0 (same month, r = 0.343) and declines monotonically across all lag steps: lag 1 (r = 0.300), lag 2 (0.283), lag 3 (0.249), lag 4 (0.212), lag 5 (0.209), and lag 6 (0.213). The absence of a peak at longer lags indicates that PDSI does directly incorporate current-month precipitation through the moisture anomaly term Z_t. However, the maximum achievable same-month correlation (0.343) is far below SPI-1 (0.76), confirming that the recursive averaging structure severely dilutes the precipitation signal. The slow decline of r across lags (remaining above 0.20 even at lag 6) reflects the long-term memory characteristic of PDSI: past precipitation months continue to exert a measurable but decreasing influence on current PDSI.

PDSI与滞后降水的中位Pearson r（图4）在同月（lag 0）最高（r = 0.343），随后单调递减：lag 1（r = 0.300）、lag 2（0.283）、lag 3（0.249）、lag 4（0.212）、lag 5（0.209）、lag 6（0.213）。无峰值出现在较长滞后处，这表明PDSI通过水分异常项Z_t确实纳入了当月降水。然而，同月可达的最大相关系数（0.343）远低于SPI-1（0.76），证实了递归平均结构对降水信号的严重稀释。r随滞后时步缓慢下降（lag 6时仍高于0.20）反映了PDSI的长期记忆特性：过去月份的降水对当期PDSI仍具有可量化但递减的影响。

### 4.6 SPI-1 vs. SPEI-1 Comparison / SPI-1与SPEI-1对比分析（Figure 5 / 图5）

The scatter plot (Figure 5, left) shows SPI-1 exceeds SPEI-1 in nearly all 124 watersheds — all points lie above the 1:1 line. The histogram of differences r(SPI-1) − r(SPEI-1) (Figure 5, right) has a median of +0.038 and is entirely right of zero, meaning only 1 watershed (0.8%) shows SPEI-1 outperforming SPI-1. This uniformity confirms that, when the objective is simply tracking same-month precipitation, SPI-1 consistently outperforms SPEI-1 across all U.S. climate regimes represented in this dataset.

散点图（图5左）显示，SPI-1在全部124个流域中几乎均优于SPEI-1，所有点均位于1:1参考线上方。差值直方图r(SPI-1) − r(SPEI-1)（图5右）中位数为+0.038，且全部分布在零值右侧——仅1个流域（0.8%）中SPEI-1优于SPI-1。这一高度一致的结果证实，当研究目标是同月降水的即时追踪时，SPI-1在本数据集覆盖的所有美国气候区内均系统性优于SPEI-1。

---

## 5. Discussion / 讨论

### 5.1 Why Does SPI-1 Correlate Best with Same-Month Rainfall? / 为何SPI-1同月相关性最高？

SPI-1's top ranking is mathematically expected. SPI-1 is a monotonic transformation of 1-month accumulated precipitation: the Gamma CDF maps precipitation to a probability in [0,1], which is then inverse-normal-transformed. Because the transformation is monotonic, the rank ordering of monthly precipitation values is preserved exactly in SPI-1. Consequently, the Spearman rank correlation ρ(SPI-1, precipitation) theoretically approaches 1.0. Pearson r is slightly below 1.0 only because the nonlinear Gamma-to-normal transformation introduces mild departures from strict proportionality.

SPI-1的最高排名在数学上是必然的。SPI-1是同月降水的单调变换：Gamma CDF将降水映射至[0,1]概率空间，再经逆正态变换得到标准正态分位数。由于变换的单调性，SPI-1完整保留了月降水值的秩次顺序。理论上，Spearman秩相关ρ(SPI-1, 降水) → 1.0；Pearson r略低于1.0，仅因非线性Gamma→正态变换引入了轻微非比例性偏差。

The implication is that SPI-1 should not be interpreted as a better drought monitor than SPEI-1, but as the most direct statistical proxy for precipitation anomalies. Whether a precipitation anomaly constitutes a hydrological drought depends on evaporative demand and antecedent conditions — considerations that SPI deliberately excludes.

因此，SPI-1不应被解读为比SPEI-1更好的干旱监测指数，而是降水异常的最直接统计代理。降水异常是否构成水文干旱取决于蒸发需求和前期条件——这些因素被SPI有意排除在外。

### 5.2 Why Does PDSI Correlate Worst with Same-Month Rainfall? / 为何PDSI同月相关性最低？

PDSI's low same-month correlation follows directly from its autoregressive structure. The recursive update equation PDSI_t = 0.897 · PDSI_{t−1} + Z_t/3 means approximately 90% of the current PDSI value is inherited from the previous month, and only ~10% reflects the current month's moisture anomaly. The time constant of τ ≈ 9.5 months means that 63% of any single-month precipitation signal is "forgotten" within approximately 9.5 months.

PDSI较低的同月相关性直接源于其自回归结构。递归更新方程PDSI_t = 0.897 · PDSI_{t−1} + Z_t/3意味着当月PDSI约90%继承自前一月，仅约10%反映当月水分异常。约9.5个月的时间常数意味着单月降水信号的63%将在约9.5个月内被"遗忘"。

As a result, PDSI at any given month is essentially a weighted average of precipitation anomalies over the preceding 6–12 months with exponentially decreasing weights. This low-pass filter behavior makes PDSI well-suited for tracking cumulative, long-duration moisture deficits but a poor indicator of instantaneous precipitation anomalies.

因此，任意时刻的PDSI本质上是前6–12个月降水异常的指数加权平均。这种低通滤波特性使PDSI非常适合追踪累积、长历时水分亏缺，但不适合作为瞬时降水异常的指示指标。

### 5.3 Why Is SPEI More Appropriate for Operational Drought Monitoring? / 为何SPEI更适合业务化干旱监测？

Despite SPI-1's higher numerical correlation with precipitation, SPEI-1 is preferable for operational drought detection because low precipitation does not always equal drought. Consider a winter month where temperatures remain below freezing: PET approaches zero, and near-zero precipitation imposes minimal water stress because evaporative demand is also near zero. SPI-1 would classify such a month as drought if precipitation falls below median, while SPEI-1 correctly reflects near-normal water balance (D = P − PET ≈ 0 − 0 ≈ 0).

尽管SPI-1与降水的数值相关性更高，但SPEI-1在业务化干旱监测中更为合适，因为低降水并不总等于干旱。以冬季气温持续低于冰点的月份为例：PET接近于零，接近零降水对景观的水分胁迫极小，因为蒸发需求同样接近零。若降水低于中位值，SPI-1会将其判定为干旱，而SPEI-1则正确反映接近正常的水量平衡状态（D = P − PET ≈ 0 − 0 ≈ 0）。

Under climate warming, rising PET intensifies drought conditions even when precipitation remains unchanged — a signal SPEI captures but SPI cannot. Vicente-Serrano et al. (2010) demonstrated that SPEI detected emerging drought trends during 1970–2008 in multiple regions that SPI failed to register.

在气候变暖背景下，即使降水不变，PET上升也会加剧实际干旱——这一信号被SPEI捕捉而SPI无法反映。Vicente-Serrano et al.（2010）证明，在1970–2008年间，SPEI在多个地区检测到了SPI未能识别的新兴干旱趋势。

### 5.4 Limitations of the Thornthwaite PET Method / PET计算方法的局限（Thornthwaite）

The Thornthwaite (1948) method requires only monthly mean temperature and latitude. However, it systematically underestimates PET in arid/windy climates and may overestimate it in cool, humid climates. Because these biases are largely systematic within each climate region, they tend to cancel in the correlation framework applied here, primarily affecting the absolute magnitude of SPEI values rather than the spatial pattern of correlation performance. Future work should evaluate sensitivity by substituting Penman-Monteith PET computed from ERA5 reanalysis data.

Thornthwaite（1948）方法仅需月均气温和纬度，但在干旱/多风气候下系统性低估PET，在寒冷湿润气候下可能高估PET。由于这些偏差在各气候区内基本为系统性偏差，在本文采用的相关性分析框架中倾向于相消，主要影响SPEI值的绝对量级而非相关性表现的空间格局。未来研究应评估用ERA5再分析数据计算Penman-Monteith PET替代后的敏感性。

### 5.5 Effect of Accumulation Scale / 时间尺度的影响（1月 vs 3月）

The decrease in correlation from 1-month to 3-month indices (SPI-1 to SPI-3: Δr ≈ −0.28, median 0.76 → 0.48; SPEI-1 to SPEI-3: Δr ≈ −0.24, median 0.69 → 0.45) reflects temporal smoothing from the longer accumulation window. The 3-month scale is not inferior in absolute terms — it is well-suited for seasonal drought monitoring, being less sensitive to single-month outliers and providing a more stable characterization of extended moisture deficits.

从1个月到3个月指数的相关性下降（SPI-1至SPI-3：Δr ≈ −0.28，中位数0.76 → 0.48；SPEI-1至SPEI-3：Δr ≈ −0.24，中位数0.69 → 0.45）反映了更长积累窗口引入的时间平滑效应。3个月尺度在绝对意义上并不逊色——它对单月异常值不敏感，能更稳定地表征持续水分亏缺，非常适合季节性干旱监测。

### 5.6 Data Limitations / 数据局限性

1. **Spatial resolution / 空间分辨率**: CMAP and Dai PDSI at 2.5° are coarse relative to many HUC4 polygons, especially in the eastern U.S. Higher-resolution alternatives (PRISM 4 km, NCLIMGRID 0.25°) would improve accuracy for small basins. / 2.5°分辨率的CMAP和Dai PDSI相对于许多HUC4多边形（尤其是东部流域）较粗糙，PRISM（4 km）或NCLIMGRID（0.25°）等高分辨率替代品将改善小流域的精度。

2. **PDSI truncation / PDSI数据截断**: PDSI ends in 2018, yielding an 18-year window versus 21 years for SPI/SPEI. / PDSI数据截至2018年，分析窗口比SPI/SPEI短约3年。

3. **Area weighting / 面积加权**: Simple arithmetic mean was used without cosine-latitude area weighting, introducing ~5–15% bias for high-latitude watersheds. / 空间聚合采用简单算术均值，未作余弦纬度面积加权，对高纬度流域引入约5–15%的误差。

4. **Calibration period overlap / 校准期重叠**: SPI/SPEI distributions were calibrated on the same period as the analysis window (2000–2020), which is appropriate for this index-comparison study but should be separated in operational forecasting applications. / SPI/SPEI分布校准期与分析窗口相同（2000–2020），适用于本指数比较研究，但在业务预报应用中应分离校准期与验证期。

---

## 6. Conclusions and Recommendations / 结论与建议

### 6.1 Main Findings / 主要发现

This study quantified the Pearson correlation between five drought indices and same-month precipitation across 124 CONUS HUC4 watersheds (2000–2020). Three principal findings emerge:

本研究量化了五种干旱指数与124个CONUS HUC4流域同月降水的Pearson相关性（2000–2020年）。三项主要发现如下：

1. **SPI-1 is the strongest concurrent precipitation proxy**, achieving median r ≈ 0.76 with near-universal statistical significance. This is theoretically expected, as SPI-1 is a monotonic transformation of 1-month precipitation.
   **SPI-1是最强的同期降水代理指标**，中位r ≈ 0.76，几乎在所有流域均达到统计显著性。这在理论上是必然结果，因为SPI-1是单月降水的单调变换。

2. **SPEI-1 is the recommended operational drought index** (median r ≈ 0.69) because it accounts for evaporative demand, avoiding false drought classifications in cold/low-PET months, and is more sensitive to climate change–driven drought intensification.
   **SPEI-1是推荐的业务化干旱监测指数**（中位r ≈ 0.69），因为它纳入了蒸发需求，避免了寒冷/低PET月份的干旱误报，且对气候变化驱动的干旱加剧更为敏感。

3. **PDSI has the lowest same-month correlation** (median r ≈ 0.34), with r declining monotonically from 0.34 at lag 0 to 0.21 at lag 6 — no peak emerges at longer lags. Its recursive averaging structure dilutes the same-month precipitation signal and imparts a long-term memory that makes it unsuitable as a concurrent precipitation proxy but valuable for tracking cumulative, multi-month moisture deficits.
   **PDSI同月相关性最低**（中位r ≈ 0.34），r从lag 0（0.34）单调递减至lag 6（0.21），无较长滞后峰值出现。其递归平均结构稀释了同月降水信号并赋予长期记忆，使其不适合作为同期降水代理，但适合追踪累积的多月水分亏缺。

### 6.2 Index Recommendations by Application / 对不同应用场景的指数建议

**Table 3 / 表3.** Recommended drought index by application context. / 不同应用场景的推荐干旱指数。

| Application / 应用场景 | Recommended Index / 推荐指数 | Rationale / 原因 |
|---|---|---|
| Real-time drought monitoring / 实时干旱监测 | **SPEI-1** | Incorporates PET; avoids false drought in cold/low-evaporation months / 纳入PET；避免寒冷/低蒸散月份误报 |
| Precipitation anomaly tracking / 降水异常监测 | SPI-1 | Direct precipitation proxy; simple computation / 直接降水代理；计算简单 |
| Agricultural & water-supply drought / 农业及供水干旱追踪 | PDSI | Captures cumulative soil moisture deficit / 捕捉累积土壤水分亏缺 |
| Seasonal drought monitoring / 季节性干旱监测 | SPI-3 / SPEI-3 | Smoothed signal; less sensitive to single-month anomalies / 平滑信号；对单月异常不敏感 |
| Climate change attribution / 气候变化归因 | SPEI | Rising PET under warming captured by SPEI / 升温驱动的PET上升被SPEI捕捉 |

### 6.3 Future Work / 后续工作展望

1. **Hydrological model calibration / 水文模型率定**: Integrate USGS stream gauge and reservoir inflow/outflow records to calibrate watershed-scale hydrological models, with SPEI-1 as a candidate lagged predictor in inflow forecasting models. / 纳入USGS流量站和水库出入库记录，以SPEI-1作为入流预报模型的候选滞后预报因子，率定流域尺度水文模型。

2. **S2S forecast integration / S2S预报集成**: Use S2S precipitation and temperature forecast products (ECMWF SEAS5, NCEP CFSv2) to compute forecast drought indices and evaluate their added predictive value relative to observed indices. / 利用S2S降水和气温预报产品（ECMWF SEAS5、NCEP CFSv2）计算预报干旱指数，评估其相对于观测指数的预测附加值。

3. **Case study watersheds / 重点流域案例分析**: Conduct depth analyses on HUC4 basins that experienced significant water-supply droughts during 2000–2020 (e.g., Colorado River headwaters, Rio Grande, Central Valley), incorporating higher-resolution data and reservoir storage records. / 对2000–2020年间经历显著供水干旱的重点流域（如科罗拉多河上游、里奥格兰德、中央谷地）开展深度案例分析，引入高分辨率数据和水库蓄水记录。

---

## References / 参考文献

Allen, R.G., Pereira, L.S., Raes, D., & Smith, M. (1998). *Crop evapotranspiration: Guidelines for computing crop water requirements* (FAO Irrigation and Drainage Paper No. 56). Food and Agriculture Organization of the United Nations.

Dai, A. (2011). Characteristics and trends in various forms of the Palmer Drought Severity Index during 1900–2008. *Journal of Geophysical Research: Atmospheres, 116*, D12115. https://doi.org/10.1029/2010JD015541

Hamlet, A.F., & Lettenmaier, D.P. (1999). Columbia River streamflow forecasting based on ENSO and PDO climate signals. *Journal of Water Resources Planning and Management, 125*(6), 333–341.

McKee, T.B., Doesken, N.J., & Kleist, J. (1993). The relationship of drought frequency and duration to time scales. *Proceedings of the 8th Conference on Applied Climatology* (pp. 179–184). American Meteorological Society.

NOAA National Centers for Environmental Information. (2023). *Billion-dollar weather and climate disasters*. National Oceanic and Atmospheric Administration.

Palmer, W.C. (1965). *Meteorological drought* (Research Paper No. 45). U.S. Weather Bureau.

Robertson, A.W., Kumar, A., Peña, M., & Vitart, F. (2015). Improving and promoting subseasonal to seasonal prediction. *Bulletin of the American Meteorological Society, 96*(3), ES49–ES53.

Thornthwaite, C.W. (1948). An approach toward a rational classification of climate. *Geographical Review, 38*(1), 55–94.

Vicente-Serrano, S.M., Beguería, S., & López-Moreno, J.I. (2010). A multiscalar drought index sensitive to global warming: The Standardized Precipitation Evapotranspiration Index. *Journal of Climate, 23*(7), 1696–1718. https://doi.org/10.1175/2009JCLI2909.1

Xie, P., & Arkin, P.A. (1997). Global precipitation: A 17-year monthly analysis based on gauge observations, satellite estimates, and numerical model outputs. *Bulletin of the American Meteorological Society, 78*(11), 2539–2558.

Yuan, X., Wood, E.F., & Ma, Z. (2015). A review on climate-model-based seasonal hydrologic forecasting: Physical understanding and system development. *WIREs Water, 2*(5), 523–536.

---

## Tables and Figures / 表格与图件

*Note: All figures are generated from the analysis notebook. Actual figure images should be inserted here.*
*注：所有图件均由分析notebook生成，实际图片应插入以下位置。*

---

**Figure 1 / 图1.** Choropleth maps of Pearson r between each of the five drought indices (SPI-1, SPEI-1, SPI-3, SPEI-3, PDSI) and same-month precipitation for ~222 CONUS HUC4 watersheds (2000–2020). Color scale: RdYlGn diverging palette; green = strong positive correlation (r > 0.8); red = near-zero or negative correlation.

五种干旱指数（SPI-1、SPEI-1、SPI-3、SPEI-3、PDSI）与约222个CONUS HUC4流域同月降水Pearson r的分级色彩地图（2000–2020年）。色带采用RdYlGn发散配色：绿色表示强正相关（r > 0.8），红色表示接近零或负相关。*(notebook Step 7.3)*

---

**Figure 2 / 图2.** Side-by-side box plots of per-HUC4 Pearson r distributions for all five drought indices. Box spans IQR (25th–75th percentile); horizontal line = median; whiskers extend to 1.5× IQR; outliers shown as individual points.

全部五种干旱指数逐HUC4 Pearson r分布的并排箱型图。箱体跨越四分位距（25th–75th百分位）；横线为中位数；须延伸至1.5倍IQR；异常值以单独点标记。*(notebook Step 7.2)*

---

**Figure 3 / 图3.** Monthly time series of precipitation and five drought index values for a representative Great Lakes HUC4 watershed (HUC4 0424, SPI-1 r = 0.76, closest to the sample median), 2000–2020. Upper panel: precipitation (mm/month). Lower panel: standardized drought index values for SPI-1, SPI-3, SPEI-1, SPEI-3, and PDSI. Dashed line marks the drought threshold (index = −0.5).

代表性大湖区HUC4流域（HUC4 0424，SPI-1 r = 0.76，最接近样本中位数）2000–2020年月降水与五种干旱指数时间序列图。上图：降水量（mm/月）；下图：SPI-1、SPI-3、SPEI-1、SPEI-3和PDSI的标准化干旱指数值。虚线为干旱阈值（指数 = −0.5）。*(generate_figures_3_4_5.py)*

---

**Figure 4 / 图4.** Median Pearson r between PDSI and lagged precipitation (lag 0–6 months) across all 124 CONUS HUC4 watersheds (2000–2018). Error bars show the 25th–75th percentile range. Lag 0 = same-month (r = 0.343); correlation declines monotonically to r ≈ 0.21 at lag 6. The darker bar highlights lag 0 as the peak.

PDSI与滞后0–6个月降水的全部124个CONUS HUC4流域中位Pearson r柱状图（2000–2018）。误差条表示各HUC4的25th–75th百分位范围。滞后0为同月（r = 0.343）；相关性单调递减至lag 6约r ≈ 0.21。深色柱标记峰值lag 0。*(generate_figures_3_4_5.py)*

---

**Figure 5 / 图5.** Left: scatter plot of r(SPI-1) versus r(SPEI-1) per HUC4 with 1:1 reference line; color indicates r(SPI-1) − r(SPEI-1). Right: histogram of differences; SPI-1 outperforms SPEI-1 in 123 of 124 watersheds (99.2%), with a median difference of +0.038.

左图：逐HUC4的r(SPI-1)与r(SPEI-1)散点图，附1:1参考线；颜色表示r(SPI-1) − r(SPEI-1)。右图：差值直方图；SPI-1在124个流域中的123个（99.2%）优于SPEI-1，中位差值为+0.038。*(generate_figures_3_4_5.py)*

---

*Abstract word count / 摘要字数: ~240 words*
*Estimated total length / 估计总篇幅: ~20 pp. double-spaced / 双倍行距约20页*
*Format / 格式: APA 7th Edition*
*Deadline / 截止日期: 2026-05-15*
