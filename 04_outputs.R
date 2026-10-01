# 导出必要的汇总结果与科学诊断图。readr 写 CSV，ggplot2/基础图形设备保存图。
# 不公开患者 ID、逐人权重或插补数据；analysis_results.rds 仅在本地使用。
as_table <- function(x, row_label = "variable") {
  d <- as.data.frame(x, stringsAsFactors = FALSE, check.names = FALSE)
  data.frame(setNames(list(rownames(d)), row_label), d, row.names = NULL, check.names = FALSE)
}
vif_output <- as_table(as.matrix(collinearity$vif))
names(vif_output)[-1] <- if (ncol(vif_output) == 2L) "VIF" else c("GVIF", "Df", "GVIF_adjusted")
public_tables <- list(
  sample_flow = bind_rows(preparation_flow, sample_flow), simulation_summary = simulation_summary,
  cleaning_summary = cleaning_summary, missingness = missingness, missing_patterns = missing_patterns,
  missing_selection = missing_selection, cc_selection = cc_selection,
  table1 = as_table(table1_output), cc_comparison = as_table(cc_comparison_output),
  collinearity = vif_output, mice_methods = mice_methods,
  mice_diagnostics = mice_postcheck, mice_mce = mice_mce,
  mice_convergence = as.data.frame(mice_convergence),
  positivity = positivity, design_summary = design_summary,
  balance_summary = balance |> group_by(design, estimand, feature) |>
    summarise(median_abs_smd = median(abs_smd), worst_abs_smd = max(abs_smd),
      worst_ecdf = max(ecdf_difference), .groups = "drop"),
  weight_diagnostics = weight_diagnostics,
  matching_quality = matching_quality, pooled_effects = pooled_effects,
  bootstrap = bootstrap_results, sensitivity_populations = sensitivity_populations
)
for (nm in names(public_tables)) readr::write_csv(public_tables[[nm]], file.path(result_dir, paste0(nm, ".csv")))
# 为可读性补充模型规格，包含预冻结的样条结点；本文件不是另一个模型配置系统。
model_specification <- data.frame(model = c("PS1", "PS2", "OR1", "OR2"),
  formula = vapply(formulas[c("ps1", "ps2", "or1", "or2")], function(f) paste(deparse(f), collapse = " "), character(1)))
readr::write_csv(model_specification, file.path(result_dir, "model_specification.csv"))
run_metadata <- data.frame(setting = c(names(config), "R_version", "causaldata_version", "completed_at_UTC",
                                      "mice_logged_events", "design_rank", "design_columns", "scaled_condition_number"),
  value = c(vapply(config, function(x) paste(x, collapse = ";"), character(1)),
            R.version.string, as.character(packageVersion("causaldata")), format(Sys.time(), tz = "UTC", usetz = TRUE),
            as.character(nrow(mice_events)), as.character(collinearity$rank), as.character(collinearity$columns),
            as.character(collinearity$scaled_condition_number)))
readr::write_csv(run_metadata, file.path(result_dir, "run_metadata.csv"))
readr::write_csv(data.frame(package = names(package_purpose),
  version = vapply(names(package_purpose), function(x) as.character(packageVersion(x)), character(1)),
  purpose = unname(package_purpose)), file.path(result_dir, "package_versions.csv"))
code_files <- c("00_main.R", "00_prepare.R", "01_missing_baseline.R", "02_ps_design.R",
                "03_effect_sensitivity.R", "functions.R", "04_outputs.R")
readr::write_csv(data.frame(file = code_files, md5 = code_md5(file.path(nhefs_dir, code_files))),
                 file.path(result_dir, "code_checksums.csv"))
input_files <- c("nhefs_rwe_dirty.csv", "nhefs_rwe_clean.csv", "nhefs_rwe_analysis_candidate.csv")
readr::write_csv(data.frame(file = input_files, md5 = code_md5(file.path(data_dir, input_files))),
                 file.path(result_dir, "input_checksums.csv"))
capture.output(sessionInfo(), file = file.path(result_dir, "sessionInfo.txt"))
saveRDS(nhefs_results, file.path(result_dir, "analysis_results.rds"))
saveRDS(imp, file.path(result_dir, "mice_object.rds")) # 本地复核轨迹/分布；.gitignore 排除。

if (config$draw_plots) {
  for (nm in names(plots)) {
    height <- if (nm == "balance") 10 else if (nm == "effects") 12 else if (nm == "mice_trace") 10 else 6
    path <- file.path(figure_dir, paste0(nm, ".png"))
    if (inherits(plots[[nm]], "ggplot")) {
      ggplot2::ggsave(path, plots[[nm]], width = 11, height = height, dpi = 160, bg = "white")
    } else {
      grDevices::png(path, width = 1760, height = height * 160, res = 160)
      tryCatch(print(plots[[nm]]), finally = grDevices::dev.off())
    }
  }
}
message("输出已保存：", result_dir, "；完成后可 quarto render report.qmd")
