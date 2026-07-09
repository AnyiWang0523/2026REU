# 报告大纲 / Report Outline

**题目 / Title:**
评估三种干旱指数与美国 HUC4 流域降水相关性的比较研究
*Comparative Assessment of Drought Index–Precipitation Correlation Across U.S. HUC4 Watersheds*

**截止日期 / Deadline:** 2026-05-15
**项目 / Program:** 2026 REU — University of Illinois Urbana-Champaign

---

## 一、引言 / 1. Introduction

**1.1 研究背景与动机**
*1.1 Background and Motivation*

- 干旱对灌溉水库运营的影响及次季节-季节（S2S）预报的重要性
  *(Drought impacts on reservoir operations; value of sub-seasonal to seasonal forecasting)*
- 近年来水库运营模型与 S2S 气象预报的显著进展
  *(Recent advances in reservoir operation modeling and S2S meteorological forecasting)*
- 在气候多样地区建立流域水文模型、理解上述进展对干旱减缓价值的重要性
  *(Importance of developing watershed-scale hydrologic models in climatologically diverse regions)*
- 干旱监测中选择合适指数的必要性
  *(Need to select appropriate drought indicators for operational use)*

**1.2 研究目标**
*1.2 Research Objectives*

- 核心问题：哪种干旱指数与同月降水的相关性最高？
  *(Core question: Which drought index most immediately tracks monthly precipitation variability?)*
- 覆盖范围：美国本土（CONUS）约 222 个 HUC4 流域单元
  *(Spatial scope: ~222 CONUS HUC4 watersheds)*
- 分析时段：2000–2020 年
  *(Analysis period: 2000–2020)*

**1.3 研究意义**
*1.3 Significance*

- 为水库入流预报模型提供输入指标选择依据
  *(Informs index selection for reservoir inflow forecasting models)*
- 为后续水文模型率定提供干旱表征框架
  *(Provides drought characterization framework for subsequent hydrological model calibration)*

---

## 二、数据与研究区域 / 2. Data and Study Area

**2.1 研究区域**
*2.1 Study Area*

- USGS WBD HUC4 流域边界，覆盖美国本土全域（HUC4 编码 01–18）
  *(USGS Watershed Boundary Dataset HUC4 units, CONUS, codes 01–18)*
- 图示：CONUS HUC4 分布地图
  *(Figure: Map of HUC4 watershed units)*

**2.2 气候数据**
*2.2 Climate Data*

| 数据集 / Dataset | 来源 / Source | 分辨率 / Resolution | 时段 / Period |
|---|---|---|---|
| 月降水 CMAP Enhanced | NOAA PSL OPeNDAP | 2.5°, 月 / monthly | 2000–2020 |
| 月气温 GHCN+CAMS | NOAA PSL OPeNDAP | 0.5°, 月 / monthly | 2000–2020 |
| PDSI (Dai scPDSI) | NOAA PSL OPeNDAP | 2.5°, 月 / monthly | 2000–2018 |

- 数据均通过 OPeNDAP 远程访问，无需下载大文件（每变量 < 2 MB）
  *(All datasets accessed remotely via OPeNDAP — no large local downloads required)*

**2.3 分析时段说明**
*2.3 Notes on Analysis Period*

- 统一时段 2000–2020，与现有水资源使用数据（WSDI）对齐
  *(Unified period 2000–2020, aligned with existing WSDI water-use records)*
- PDSI 数据截至 2018 年，对应分析窗口缩短至 2000–2018
  *(PDSI dataset ends in 2018; PDSI correlation window shortened accordingly)*

---

## 三、方法 / 3. Methods

**3.1 数据获取与空间聚合**
*3.1 Data Retrieval and Spatial Aggregation* (`src/aggregate.py`)

- 通过 OPeNDAP 远程获取 NOAA PSL 月尺度数据（选用 pydap 引擎避免 Windows 平台 libcurl 兼容问题）
  *(Remote OPeNDAP access using pydap engine — avoids libcurl issues on Windows)*
- 经纬度坐标标准化：0–360° 经度转换为 −180–180°；纬度由南→北排序
  *(Coordinate normalization: 0–360° longitude → −180–180°; latitude sorted ascending)*
- CMAP 单位转换：mm/day × 月天数 → mm/month
  *(Unit conversion: mm/day × days_in_month → mm/month)*
- regionmask 建立 HUC4 掩膜，逐流域空间均值聚合（lat × lon → n_HUC4）
  *(Spatial aggregation: regionmask polygon mask + groupby mean → per-HUC4 time series)*

**3.2 干旱指数计算**
*3.2 Drought Index Computation* (`src/compute_indices.py`)

**SPI — 标准化降水指数（McKee et al. 1993）**
*SPI — Standardized Precipitation Index*

1. 滑动窗口累加降水（1 个月或 3 个月）
   *(Rolling accumulation of precipitation over 1- or 3-month window)*
2. 对校准期（2000–2020）逐历月（Jan–Dec）拟合 Gamma 分布（含零降水混合分布）
   *(Gamma distribution fit per calendar month, mixed with zero-precipitation probability)*
3. Probit 变换：SPI = Φ⁻¹(H(P_acc))，标准化为 N(0,1)
   *(Probit transform to standard normal)*

**SPEI — 标准化降水蒸散指数（Vicente-Serrano et al. 2010）**
*SPEI — Standardized Precipitation-Evapotranspiration Index*

1. 月水分平衡：D = P − PET（结合降水与潜在蒸散）
   *(Monthly water balance: D = P − PET, combining precipitation and evaporative demand)*
2. 滑动窗口累加 D（1 个月或 3 个月）
   *(Rolling accumulation of D over 1- or 3-month window)*
3. 拟合 Log-Logistic 分布（由 Pearson III 近似），再做 Probit 变换
   *(Log-logistic distribution fit via Pearson III approximation, then probit transform)*
4. PET 采用 Thornthwaite（1948）法：仅需月均气温与纬度（适合数据有限场景）
   *(PET computed by Thornthwaite method using temperature and latitude only)*

**PDSI — 帕尔默干旱烈度指数（Palmer 1965, Dai 自校准版）**
*PDSI — Palmer Drought Severity Index*

1. 双层土壤水量平衡（表层 + 深层）
   *(Two-layer soil water balance: surface + deep layer)*
2. 递归平滑：PDSI_t = 0.897 · PDSI_{t-1} + Z/3
   *(Recursive smoothing: inherent multi-month memory)*
3. 时间常数约 10 个月，具有 3–6 个月的滞后记忆效应
   *(Time constant ~10 months; responds to cumulative past conditions)*
4. 直接使用 Dai scPDSI 产品（NOAA PSL），不重新计算
   *(Dai scPDSI product used directly from NOAA PSL)*

**3.3 相关性分析**
*3.3 Correlation Analysis* (`src/correlation.py`)

- 逐 HUC4 计算各指数与同月降水的 Pearson r
  *(Per-HUC4 Pearson r between each drought index and same-month rainfall)*
- 双侧 p 值检验（显著性阈值 α = 0.05）
  *(Two-sided significance test, α = 0.05)*
- PDSI 滞后分析：将降水分别提前 0–6 个月，计算与 PDSI 的中位相关系数
  *(PDSI lag analysis: rainfall shifted 0–6 months back, median r computed per lag)*

---

## 四、结果 / 4. Results

**4.1 各指数与同月降水相关性汇总**
*4.1 Summary Statistics of Index–Rainfall Correlation*

- 表格：各指数中位数 r、均值 r、标准差、显著流域比例（p < 0.05）
  *(Table: median r, mean r, std r, fraction of significant HUC4s per index)*
- 排序结果（由高到低）：SPI-1 > SPEI-1 > SPI-3 > SPEI-3 > PDSI
  *(Ranking, highest to lowest: SPI-1 > SPEI-1 > SPI-3 > SPEI-3 > PDSI)*

**4.2 空间分布图**
*4.2 Spatial Maps*

- **图 1**：五个指数 Pearson r 的 HUC4 空间分布图（RdYlGn 色带；绿色=强正相关，红色=负相关）
  *(Figure 1: Choropleth maps of Pearson r per HUC4 for all five indices, RdYlGn colormap)*

**4.3 各 HUC4 相关系数分布箱型图**
*4.3 Box Plot of Pearson r Across HUC4s*

- **图 2**：五个指数的箱型图并排对比，直观展示中位数与离散程度
  *(Figure 2: Side-by-side box plots comparing r distribution across all five indices)*

**4.4 单个流域时间序列示例**
*4.4 Sample Time Series (One HUC4)*

- **图 3**：选定流域的月降水与各指数时间序列叠加图（2000–2020）
  *(Figure 3: Monthly time series of rainfall and all drought indices for a selected watershed)*

**4.5 PDSI 滞后分析**
*4.5 PDSI Lag Analysis*

- **图 4**：PDSI 与滞后 0–6 月降水的中位 Pearson r 柱状图
  *(Figure 4: Bar chart of median r between PDSI and rainfall shifted 0–6 months)*
- 结果：峰值相关出现在约 t−2 至 t−3 月，印证 PDSI 的滞后记忆效应
  *(Finding: peak correlation at ~2–3 month lag, confirming PDSI's memory effect)*

**4.6 SPI 与 SPEI 的差异分析**
*4.6 SPI vs. SPEI Comparison*

- **图 5**：SPI-1 vs SPEI-1 散点图（逐 HUC4）及 r(SPI-1) − r(SPEI-1) 差值直方图
  *(Figure 5: Scatter plot of r(SPI-1) vs r(SPEI-1) per HUC4, and histogram of differences)*

---

## 五、讨论 / 5. Discussion

**5.1 为何 SPI-1 同月相关性最高？**
*5.1 Why Does SPI-1 Correlate Best with Same-Month Rainfall?*

- SPI-1 是同月降水的单调变换（Gamma CDF → 正态分位数），秩次顺序与降水完全一致
  *(SPI-1 is a monotonic transform of same-month precipitation — rank order is preserved)*
- 理论上 Spearman ρ(SPI-1, rainfall) → 1.0；Pearson r 略低于 1 仅因变换非线性
  *(Theoretically Spearman ρ → 1.0; Pearson r slightly below 1 due to nonlinear transform)*

**5.2 为何 PDSI 同月相关性最低？**
*5.2 Why Does PDSI Correlate Worst with Same-Month Rainfall?*

- 递归系数 0.897 产生约 10 个月的时间常数（低通滤波效应）
  *(Recursive coefficient 0.897 creates ~10-month time constant — acts as a low-pass filter)*
- 当月降水仅贡献约 10% 的 PDSI 变化，其余来自历史积累
  *(Same-month precipitation contributes only ~10% of current PDSI; the rest is inherited history)*
- 因此 PDSI 不适合作为同月降水的代理指标，但适合追踪累积水分亏缺
  *(PDSI is a poor same-month proxy but ideal for tracking cumulative moisture deficit)*

**5.3 为何 SPEI 更适合实际干旱监测？**
*5.3 Why Is SPEI More Appropriate for Real-Time Drought Monitoring?*

- 低降水 ≠ 干旱：当气温低、PET 极小时，0 mm 降水并不造成实际水分亏缺
  *(Low precipitation ≠ drought: when PET is near zero in cold months, no-rain ≠ drought)*
- 典型场景：若某月前后两个月均有降水，而该月无降水但气温低（冬/春季），SPI-1 会误报干旱；SPEI-1 通过水分平衡 D = P − PET 正确识别该月水分充足
  *(Example: a no-rain month sandwiched between wet months in winter is not a drought — SPEI-1 correctly reflects near-normal water balance while SPI-1 misclassifies it as drought)*
- 在气候变暖背景下，PET 上升使实际干旱加剧，SPEI 能捕捉这一信号而 SPI 不能
  *(Under warming climate, rising PET intensifies drought — SPEI captures this, SPI cannot)*

**5.4 PET 计算方法的局限（Thornthwaite）**
*5.4 Limitations of the Thornthwaite PET Method*

- 仅依赖温度，在干旱/半干旱地区系统性低估蒸散，在湿润地区高估
  *(Temperature-only method underestimates PET in arid/windy climates, overestimates in humid ones)*
- 对相关分析影响有限（系统偏差在 HUC4 间相消），但影响绝对量级
  *(Systematic bias largely cancels in correlation analysis but affects absolute PET magnitudes)*
- 更精确的替代方案：Penman-Monteith（需辐射、风速数据）
  *(More accurate alternative: Penman-Monteith, requiring radiation and wind data)*

**5.5 时间尺度的影响（1 月 vs 3 月）**
*5.5 Effect of Accumulation Scale (1-month vs 3-month)*

- 3 个月指数引入平滑滞后，导致与同月降水相关性下降
  *(3-month accumulation introduces smoothing lag, reducing same-month correlation)*
- 3 个月尺度更适合季节性干旱监测（信号更稳定、噪声更低）
  *(3-month scale better suited for seasonal drought monitoring with reduced noise)*

**5.6 数据局限性**
*5.6 Data Limitations*

- CMAP / PDSI 分辨率为 2.5°，对面积较小的 HUC4 流域代表性有限
  *(2.5° resolution limits accuracy for small HUC4 watersheds)*
- PDSI 数据截至 2018 年，分析窗口比其他指数短约 2 年
  *(PDSI ends in 2018, yielding a shorter analysis window than other indices)*
- 空间聚合采用简单算术均值，未作余弦面积加权（对高纬度流域存在约 15% 的误差）
  *(Simple arithmetic mean used in spatial aggregation; no cosine area weighting applied)*

---

## 六、结论与建议 / 6. Conclusions and Recommendations

**6.1 主要发现**
*6.1 Main Findings*

- SPI-1 与同月降水的 Pearson r 中位数最高（约 0.76），是降水异常的最直接统计代理
  *(SPI-1 achieves highest median Pearson r ~0.76 — most direct statistical proxy for precipitation)*
- SPEI-1 略低（约 0.69），但综合了蒸散需求，对气候变化更敏感
  *(SPEI-1 slightly lower ~0.69 but incorporates evaporative demand; more climate-change-sensitive)*
- PDSI 相关性最低（约 0.34），且峰值相关滞后约 2–3 个月，反映其低通滤波特性
  *(PDSI lowest ~0.34, with peak correlation at 2–3 month lag, reflecting low-pass filter behavior)*

**6.2 对不同应用场景的指数建议**
*6.2 Index Recommendations by Application*

| 应用场景 / Use Case | 推荐指数 / Recommended Index | 原因 / Reason |
|---|---|---|
| 实时干旱监测 / Real-time drought monitoring | **SPEI-1** | 纳入 PET，低降水月不一定是干旱；避免冬季/低温月误报 / Incorporates PET; low-rain months with low PET are not droughts |
| 纯降水异常监测 / Precipitation anomaly tracking | SPI-1 | 与降水直接对应，计算简单 / Direct precipitation proxy, simple to compute |
| 农业干旱 / 供水危机追踪 / Agricultural & water-supply drought | PDSI | 捕捉累积土壤水分亏缺，滞后反映长期干旱 / Captures cumulative soil moisture deficit |
| 季节性干旱监测 / Seasonal drought monitoring | SPI-3 / SPEI-3 | 平滑短期波动，更稳定的季节信号 / Smoothed signal, less sensitive to single-month anomalies |
| 气候变化归因 / 未来情景 / Climate change attribution & projections | SPEI | PET 随升温上升，SPEI 能捕捉增温驱动的干旱 / Rising PET under warming captured by SPEI |

**6.3 后续工作展望**
*6.3 Future Work*

- 引入灌溉水库出入库径流记录，建立并率定流域水文模型
  *(Incorporate reservoir inflow/outflow records for hydrological model development and calibration)*
- 接入 S2S 气象预报数据，评估预报干旱指数相对于观测指数的预测附加值
  *(Integrate S2S meteorological forecasts to assess their added value for drought prediction)*
- 对近期经历供水干旱的重点流域开展案例深度分析
  *(Conduct case studies on watersheds that experienced recent water-supply droughts)*

---

## 七、参考文献 / 7. References

- McKee, T.B., Doesken, N.J., & Kleist, J. (1993). The relationship of drought frequency and duration to time scales. *Proc. 8th Conference on Applied Climatology*, AMS.
- Palmer, W.C. (1965). *Meteorological Drought*. U.S. Weather Bureau Research Paper No. 45.
- Vicente-Serrano, S.M., Beguería, S., & López-Moreno, J.I. (2010). A multiscalar drought index sensitive to global warming: The Standardized Precipitation Evapotranspiration Index. *Journal of Climate, 23*(7), 1696–1718.
- Dai, A. (2011). Characteristics and trends in various forms of the Palmer Drought Severity Index during 1900–2008. *Journal of Geophysical Research: Atmospheres, 116*, D12115.
- Xie, P., & Arkin, P.A. (1997). Global precipitation: A 17-year monthly analysis based on gauge observations, satellite estimates, and numerical model outputs. *Bulletin of the American Meteorological Society, 78*(11), 2539–2558.
- Thornthwaite, C.W. (1948). An approach toward a rational classification of climate. *Geographical Review, 38*(1), 55–94.

---

## 附图清单 / List of Figures

| 图号 / Fig. | 内容 / Content | 对应代码 / Code Location |
|---|---|---|
| 图 1 / Fig. 1 | 各指数 HUC4 Pearson r 空间分布图（5 幅）/ Choropleth maps of r per index | notebook Step 7.3 |
| 图 2 / Fig. 2 | 各指数 r 分布箱型图 / Box plots of r across HUC4s | notebook Step 7.2 |
| 图 3 / Fig. 3 | 示例流域月时间序列图 / Sample HUC4 monthly time series | notebook Step 6.1 |
| 图 4 / Fig. 4 | PDSI 滞后相关分析柱状图 / PDSI lag correlation bar chart | notebook Step 6.2 |
| 图 5 / Fig. 5 | SPI-1 vs SPEI-1 散点图与差值直方图 / SPI-1 vs SPEI-1 scatter and diff histogram | notebook Step 8.1 |
