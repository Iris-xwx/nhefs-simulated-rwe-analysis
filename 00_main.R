# NHEFS 戒烟与长期体重变化：主体统计及因果分析
# 在项目根目录或 nhefs 目录运行 source("nhefs/00_main.R", encoding="UTF-8")
# （若已在 nhefs 目录，改为 source("00_main.R", encoding="UTF-8")）。
# 输入：causaldata::nhefs；本目录可独立完成模拟降级、确定性清洗与分析。
# 本版保留图中 16 项主体方法，默认导出必要结果；Quarto 只读取输出，不重复拟合。
# 正式默认值：20 次插补、20 次迭代、全部插补分别匹配、500 次 bootstrap。

package_purpose <- c(
  causaldata = "提供公开教学数据 nhefs；不把模拟 EHR/Registry 标签称为真实来源",
  dplyr = "数据整理、分组汇总、合并各次插补结果",
  readr = "读取当前项目的清洗 CSV；可选导出核心结果表",
  ggplot2 = "缺失图、ECDF、PS 分布、Love plot 和效应森林图",
  mice = "PMM/Logistic/多分类 MICE、插补诊断、Rubin 合并",
  tableone = "Table 1、完整案例保留/删除人群比较、基线 SMD",
  car = "设计矩阵多重共线性诊断中的 VIF/GVIF",
  splines = "固定结点的自然样条（R 自带推荐包）",
  sandwich = "异方差稳健 HC0 与匹配子类聚类稳健方差",
  MatchIt = "完全匹配、最近邻、卡尺和遗传匹配的统一接口",
  optmatch = "MatchIt 完全匹配所需的优化后端",
  Matching = "MatchIt 遗传匹配所需的匹配与优化后端",
  rgenoud = "遗传匹配搜索协变量距离权重的优化后端"
)
missing_packages <- names(package_purpose)[!vapply(
  names(package_purpose), requireNamespace, logical(1), quietly = TRUE
)]
if (length(missing_packages)) stop(
  "请先安装缺少的包：", paste(missing_packages, collapse = ", ")
)
suppressPackageStartupMessages(library(dplyr)) # 其余包使用 包名::函数，便于追溯用途

nhefs_dir <- if (file.exists("00_main.R") && file.exists("functions.R")) "." else "nhefs"
nhefs_dir <- normalizePath(nhefs_dir, winslash = "/", mustWork = TRUE)
data_dir <- file.path(nhefs_dir, "data")
result_dir <- file.path(nhefs_dir, "results")
figure_dir <- file.path(result_dir, "figures")
dir.create(data_dir, showWarnings = FALSE)
dir.create(figure_dir, showWarnings = FALSE, recursive = TRUE)

# 参数可在运行入口前用 options(nhefs.m=..., ...) 覆盖，便于小规模验证。
config <- list(
  m = getOption("nhefs.m", 20L), maxit = getOption("nhefs.maxit", 20L),
  seed = 20260928L, representative = 13L,
  run_matching = getOption("nhefs.matching", TRUE),
  bootstrap = getOption("nhefs.bootstrap", 500L), bootstrap_seed = 94113L,
  draw_plots = getOption("nhefs.plots", TRUE),
  display_plots = getOption("nhefs.display", interactive()),
  export_results = getOption("nhefs.export", TRUE),
  rebuild = getOption("nhefs.rebuild", FALSE), preparation_seed = 20260926L,
  smd_limit = 0.10, caliper = 0.20,
  genetic_pop = getOption("nhefs.genetic_pop", 50L),
  genetic_generations = getOption("nhefs.genetic_generations", 10L),
  genetic_wait = 4L, ps_trim = c(0.05, 0.95), truncation = c(0.01, 0.99)
)

stopifnot(config$m >= 2L, config$maxit >= 3L, config$bootstrap >= 0L)
source(file.path(nhefs_dir, "functions.R"), encoding = "UTF-8", local = TRUE)
source(file.path(nhefs_dir, "00_prepare.R"), encoding = "UTF-8", local = TRUE)
source(file.path(nhefs_dir, "01_missing_baseline.R"), encoding = "UTF-8", local = TRUE)
source(file.path(nhefs_dir, "02_ps_design.R"), encoding = "UTF-8", local = TRUE)
source(file.path(nhefs_dir, "03_effect_sensitivity.R"), encoding = "UTF-8", local = TRUE)

# 核心交付对象：在 RStudio Environment 中查看 nhefs_results。
# design_objects 保留每次插补的匹配权重与 subclass，避免堆叠数据后错误估计。
nhefs_results <- list(
  simulation = simulation_summary, cleaning = cleaning_summary,
  sample_flow = sample_flow, missingness = missingness, missing_patterns = missing_patterns,
  missing_selection = missing_selection, cc_selection = cc_selection,
  table1 = table1_output, collinearity = collinearity,
  mice_methods = mice_methods, mice_events = mice_events,
  mice_convergence = mice_convergence, mice_postcheck = mice_postcheck,
  mice_mce = mice_mce, positivity = positivity,
  balance = balance, design_summary = design_summary, weight_diagnostics = weight_diagnostics,
  effect_by_imputation = effect_by_imputation, pooled_effects = pooled_effects,
  bootstrap = bootstrap_results, influence = influence_diagnostics,
  sensitivity_populations = sensitivity_populations, matching_status = matching_status
)
print(pooled_effects[, c("method", "estimand", "estimate", "se", "conf_low", "conf_high")])
if (config$export_results) {
  source(file.path(nhefs_dir, "04_outputs.R"), encoding = "UTF-8", local = TRUE)
}
