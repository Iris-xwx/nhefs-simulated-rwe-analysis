# NHEFS: Simulated RWE Data Quality & Causal Analysis
# NHEFS：模拟 RWE 数据质量与因果分析

**Author / 作者:** Xuchi Cheng · [Iris-xwx](https://github.com/Iris-xwx)  
**Tools / 工具:** R · Quarto · MICE · propensity scores · causal inference  
**License / 许可:** [MIT](LICENSE)

[Read the report / 在线阅读报告](https://iris-xwx.github.io/nhefs-simulated-rwe-analysis/) ·
[中文说明](#中文说明) · [Reproduce / 复现](#reproduce--复现)

## Overview

A reproducible personal portfolio project examining smoking cessation and
1971–1982 weight change using public NHEFS observational teaching data.
It introduces controlled, RWE-style data quality problems and follows the
analysis from deterministic cleaning and missing-data handling through
propensity-score design, effect estimation and sensitivity analysis.

**This is a methods demonstration using public teaching data with simulated
quality degradation. It is not a newly collected EHR/registry study, a complete
replication of textbook results, or clinical evidence for treatment decisions.**

## What this project demonstrates

- **Data quality:** fixed-seed heterogeneous coding, sentinel values, mixed
  kg/lb units, logical conflicts and duplicate extracts; deterministic cleaning,
  quarantine rules and an auditable sample flow.
- **Missing data:** missingness and selection diagnostics, complete-case analysis,
  20 MICE imputations with 20 iterations, convergence checks and Rubin pooling.
- **Causal design:** logistic propensity scores with fixed-knot natural splines;
  positivity, common support, SMD/ECDF balance and effective sample size;
  IPTW, stabilized weights, ATT, overlap weights, trimming and truncation.
- **Matching:** full, nearest-neighbor, caliper and genetic matching in each
  imputation, with a prespecified balance gate before effect reporting.
- **Estimation:** outcome standardization, AIPW, robust/influence-function variance,
  refitted bootstrap estimates and sensitivity analyses with explicit estimands.
- **Reproducibility:** modular R scripts, a version lockfile, source/input checksums,
  recorded parameters, package versions, automated consistency checks and a
  self-contained Quarto report.

## Selected outputs

The committed run has 1,629 unique people, 1,562 analysis candidates and 1,329
complete cases. All main multiple-imputation estimates pool 20 imputations.

| Method | Target | Weight-change difference, kg | 95% CI |
|---|---|---:|---|
| IPTW | Candidate-population ATE | 3.50 | 2.48 to 4.53 |
| AIPW | Candidate-population ATE | 3.43 | 2.49 to 4.36 |

Positive estimates indicate greater weight gain in the cessation group.
These are estimates for this simulated-quality teaching analysis under the
specified identification and modeling assumptions, not unconditional clinical
causal conclusions.

Stabilized IPTW passes the `|SMD| < 0.10` balance gate in all imputations.
None of the four matching designs passes that gate in every imputation, so
only their diagnostics are shown, without formal matching effect estimates.
Different target populations are reported separately.

![Effect estimates by method and target population](results/figures/effects.png)

- [Full report (Chinese)](https://iris-xwx.github.io/nhefs-simulated-rwe-analysis/)
- [Report source](report.qmd) · [Downloadable HTML](report.html)
- [Pooled estimates](results/pooled_effects.csv)
- [Design and balance summary](results/design_summary.csv)
- [Figures](results/figures/) · [Run parameters](results/run_metadata.csv)

## 中文说明

这是一个用于展示 R 数据分析与因果推断能力的可复现个人项目。研究问题是戒烟与
1971—1982 年体重变化的调整后差异。基于公开 NHEFS 教学数据，人为加入 RWE
风格的数据质量问题，再完成清洗、缺失处理、倾向评分设计、效应估计和敏感性分析。

**项目是公开观察性教学数据上的方法演示，不是真实 EHR/Registry 研究，也不是
原始教材结果的完全复刻。人工来源标签、缺失机制及数据损坏均公开说明。**

主要展示能力：

- 数据质量处理：编码异构、哨兵值、kg/lb 混用、逻辑冲突、重复提取及可追溯清洗。
- 缺失值处理：缺失模式与选择诊断、完整案例、MICE 插补、收敛检查与 Rubin 合并。
- 因果设计：Logistic PS 与固定结点自然样条、正值性、共同支持、SMD/ECDF、
  有效样本量、IPTW/稳定权重/ATT/overlap、修剪与权重截尾。
- 匹配与估计：四种匹配、预先指定的平衡门槛、结局标准化、AIPW、稳健方差及
  重新拟合模型的 bootstrap；区分不同目标人群与 estimand。
- 工程复现：模块化脚本、版本锁定、源码与输入指纹、参数记录、结果检查与 Quarto 报告。

当前正式结果包含 1,629 名去重后的个体、1,562 名分析候选者及 1,329 名完整案例。
20 次插补合并后，IPTW 的 ATE 为 3.50 kg（95% CI 2.48–4.53），AIPW 为
3.43 kg（95% CI 2.49–4.36）。正值表示戒烟组相对体重增加更多。
稳定 IPTW 在所有插补中达到 `|SMD| < 0.10`；四种匹配方案均未在所有插补中
通过门槛，因此只展示诊断。上述估计不能解释为无条件成立的临床因果结论。

## Project structure / 项目结构

| File / 文件 | Purpose / 用途 |
|---|---|
| `00_main.R` | Analysis entry, parameters and package roles / 分析入口与参数 |
| `00_prepare.R` | Teaching inputs, simulated damage, cleaning / 数据与清洗 |
| `01_missing_baseline.R` | Missingness, MICE, baseline summaries / 缺失与基线 |
| `02_ps_design.R` | PS, weights, matching, diagnostics / 倾向评分设计 |
| `03_effect_sensitivity.R` | Estimation, pooling, sensitivity / 估计与敏感性分析 |
| `functions.R` | Shared statistical helpers / 共用统计函数 |
| `04_outputs.R` | Tables, figures and run records / 输出与环境记录 |
| `report.qmd`, `_quarto.yml`, `report.css` | Report source and presentation / 报告与样式 |
| `check_project.R` | Read-only consistency checks / 只读一致性检查 |
| `setup_dependencies.R`, `renv.lock` | Dependencies and verified versions / 依赖与版本 |
| `data/README.md` | Data provenance and local regeneration / 数据出处与再生成 |
| `results/` | Aggregate results and figures / 汇总结果与图 |
| `docs/index.html` | Published report copy for GitHub Pages / 网页发布副本 |

## Reproduce / 复现

Verified environment: **R 4.4.1 and Quarto 1.10.18**. Quarto is a separate
application. Clone/download the repository and work from its root directory.

已验证环境为 **R 4.4.1、Quarto 1.10.18**。Quarto 应用需单独安装。
克隆或下载仓库后，从仓库根目录运行。

### 1. Dependencies / 安装依赖

For a convenient setup, install missing direct dependencies:
便捷方式安装缺少的直接依赖，不锁定版本：

```r
source("setup_dependencies.R", encoding = "UTF-8")
```

For the recorded package versions, use `renv` instead:
按已记录版本复现可改用：

```r
install.packages("renv")
renv::restore(project = ".", lockfile = "renv.lock", prompt = FALSE)
# Restart R, then explicitly load the project library.
# 重启 R 后显式载入项目库。
renv::load(project = ".")
```

### 2. Analysis / 运行分析

```r
source("00_main.R", encoding = "UTF-8")
```

Defaults: `m=20`, 20 MICE iterations, four matching designs in all imputations,
and 500 bootstrap replicates. Genetic matching can take substantial time.
Person-level inputs are regenerated locally from `causaldata::nhefs`.
An existing dirty extract is reused unless `options(nhefs.rebuild=TRUE)` is set.

默认进行 20 次插补、每次 20 次 MICE 迭代、所有插补中的四种匹配及 500 次
bootstrap；遗传匹配较耗时。逐人输入由公开教学数据在本地生成。
已有 dirty CSV 默认不覆盖；需要重建时先设置 `options(nhefs.rebuild=TRUE)`。

### 3. Report and checks / 渲染与检查

```sh
quarto render report.qmd
Rscript --vanilla check_project.R
```

If R is not on your terminal path, run
`source("check_project.R", encoding="UTF-8")` in RStudio. On Windows, use a
UTF-8-capable R locale. The report reads committed aggregate outputs and does
not rerun models. Changing report text only requires rendering; changing
analysis code, inputs or parameters requires a fresh analysis.

如果 R 不在终端搜索路径中，可在 RStudio 执行上述 `source()`；Windows 下应
使用支持 UTF-8 的 R 区域设置。报告只读取汇总输出，不重跑分析。
只改文字可直接渲染；修改分析代码、输入或参数后应重新分析。

Small smoke runs may use:
小规模代码检查可临时设置：

```r
options(nhefs.m=2, nhefs.maxit=3, nhefs.matching=FALSE, nhefs.bootstrap=0)
```

These settings are not valid for the formal report. Start a fresh R session
before the formal run. Missing required outputs or changed source/input
checksums prevent formal report generation.

这些参数不用于正式报告；正式运行前启动新的 R 会话。
缺少必要输出或源码/输入指纹不一致时，不能据旧结果生成正式报告。

## Website updates / 网页更新

GitHub Pages serves the pre-rendered report from `main:/docs`. After rendering
and checking, update the published copy in R and commit/push the changes:

GitHub Pages 从 `main` 分支的 `/docs` 发布已渲染报告。完成渲染和检查后，
在 R 中同步副本，再提交并推送：

```r
stopifnot(file.copy("report.html", "docs/index.html", overwrite = TRUE))
```

The website publishes the existing report; it does not run the statistical
pipeline on GitHub. The full report is in Chinese; this README is bilingual.

网页发布现有报告，不在 GitHub 上执行统计分析。完整报告为中文，README 为中英双语。

## Interpretation and scope / 解释边界

Design uses baseline information other than the outcome; imputation may use
the observed outcome. Matching effects require the balance gate to pass in
every imputation. ATT, overlap, complete-case and restricted populations are
not presented as interchangeable replications of the candidate-population ATE.

Weighted sandwich variance conditions on fitted weights; outcome-regression
delta variance conditions on the empirical covariate distribution; same-sample
AIPW influence-function variance is approximate. Bootstrap resamples the fixed
13th representative imputation and refits PS/outcome models; it does not repeat
the full MICE pipeline and is not a full MI uncertainty interval.

Identification depends on consistency, conditional exchangeability and
positivity. The project does not fully address outcome-missingness MAR/MNAR,
precise cessation-time alignment, survey weights/complex sampling, cross-fitting
or external generalizability. Fixed seeds do not guarantee bitwise equality
across R/package versions.

设计使用结局以外的基线信息，插补可以使用已观测结局。不同目标人群单独报告。
权重模型与结局模型方差存在上述条件性与近似；bootstrap 只重抽样固定的第 13 次
代表插补并重拟合模型，不重做完整 MICE，不能当作 MI 全流程的不确定性区间。
识别依赖一致性、条件交换性与正值性。未开展完整结局缺失 MAR/MNAR 分析、
精确戒烟时点对齐、复杂抽样推断、交叉拟合或外部可推广性评估。
固定种子不保证跨版本逐位完全一致。

## Data, license and references / 数据、许可与参考资料

Original code and documentation: [MIT License](LICENSE), © 2026 Xuchi Cheng.
See [license and attribution notes](LICENSE_NOTES.md) for upstream materials.
Person-level CSVs, local RDS files, R session state, caches and credentials are
excluded by `.gitignore`; aggregate tables, figures and reports are committed.
Browser drag-and-drop uploads do not apply `.gitignore` automatically.

原创代码与文档采用 MIT 许可；上游数据与第三方材料见许可说明。
逐人 CSV、本地 RDS、R 会话状态、缓存与凭证不上传，保留汇总结果、图和报告。
浏览器拖拽上传不会自动应用 `.gitignore`。

- Hernán MA, Robins JM. *Causal Inference: What If*.
  [Book and data / 教材与数据](https://miguelhernan.org/whatifbook).
- [causaldata](https://github.com/NickCH-K/causaldata): public teaching-data input.
- [mice](https://amices.org/mice/) · [MatchIt](https://kosukeimai.github.io/MatchIt/)
  · [sandwich](https://sandwich.r-forge.r-project.org/).
- [Quarto](https://quarto.org/docs/projects/code-execution.html)
  · [renv](https://rstudio.github.io/renv/reference/snapshot.html).
