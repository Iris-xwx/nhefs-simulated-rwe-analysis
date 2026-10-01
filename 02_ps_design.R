# 6-11：PS Logistic、自然样条、正值性/共同支持、权重、修剪/截尾、四种匹配。
# 本阶段模型和设计只使用基线与暴露；结局不能参与 PS 调参或设计选择。
message("[2/3] 倾向评分与设计诊断")
formulas <- make_formulas(candidate)
design_objects <- vector("list", length(completed))
balance_parts <- weight_parts <- positivity_parts <- matching_parts <- list()
design_estimands <- c(UNADJUSTED = "Descriptive", IPTW = "ATE", STABILIZED_IPTW = "ATE",
  ATT = "ATT", OVERLAP = "Overlap", TRUNCATED_IPTW = "Modified-weight ATE",
  FULL = "ATE", NEAREST = "ATT", CALIPER = "ATT", GENETIC = "ATT")

for (i in seq_along(completed)) {
  message("  插补 ", i, "/", length(completed))
  d <- completed[[i]]
  pp <- fit_ps(d, formulas$ps2)
  p <- pp$ps; lp <- qlogis(p)
  support <- support_limits(d$z, pp$raw_ps)
  bins <- cut(p, unique(quantile(p, seq(0, 1, .1))), include.lowest = TRUE)
  decile_counts <- table(bins, factor(d$z, levels = 0:1))
  clinical_counts <- table(cut(d$age, c(-Inf, 39, 49, 59, Inf)),
                          cut(d$smokeintensity, c(-Inf, 15, 25, Inf)),
                          factor(d$z, levels = 0:1))
  positivity_parts[[i]] <- data.frame(imputation = i, converged = pp$fit$converged,
    ps_min = min(pp$raw_ps), ps_max = max(pp$raw_ps),
    support_lower = support[1], support_upper = support[2],
    outside_support = sum(p < support[1] | p > support[2]),
    ps_below_01 = sum(p < .01), ps_below_05 = sum(p < .05),
    ps_above_95 = sum(p > .95), ps_above_99 = sum(p > .99),
    empty_decile_cells = sum(decile_counts == 0), sparse_clinical_cells = sum(clinical_counts < 5))
  weights <- make_weights(d$z, p, config$truncation)
  matches <- list()
  if (config$run_matching) {
    match_args <- list(formula = formulas$ps1, data = d, distance = lp)
    specs <- list(
      FULL = list(method = "full", estimand = "ATE"),
      NEAREST = list(method = "nearest", estimand = "ATT", ratio = 1,
                     replace = FALSE, m.order = "largest"),
      CALIPER = list(method = "nearest", estimand = "ATT", ratio = 1,
                     replace = FALSE, m.order = "largest", caliper = config$caliper, std.caliper = TRUE),
      GENETIC = list(method = "genetic", estimand = "ATT", ratio = 1, replace = FALSE,
        m.order = "largest", pop.size = config$genetic_pop,
        max.generations = config$genetic_generations, wait.generations = config$genetic_wait,
        fit.func = "qqmax.max", include.obj = TRUE)
    )
    for (nm in names(specs)) {
      set.seed(73100L + i) # 遗传匹配在每个插补中独立运行，固定随机种子。
      m <- do.call(MatchIt::matchit, c(match_args, specs[[nm]]))
      matches[[nm]] <- m
      weights[[nm]] <- as.numeric(m$weights[rownames(d)])
      if (anyNA(weights[[nm]])) stop("匹配权重与患者行名无法对应")
      # 真实卡尺 = 每次插补 logit(PS) SD 的 0.2 倍；显式检查距离是否违反卡尺。
      distances <- numeric()
      if (!is.null(m$match.matrix)) {
        mm <- m$match.matrix; score <- setNames(lp, rownames(d))
        for (j in seq_len(nrow(mm))) {
          controls <- mm[j, !is.na(mm[j, ]), drop = TRUE]
          distances <- c(distances, abs(score[rownames(mm)[j]] - score[controls]))
        }
      }
      violations <- if (nm == "CALIPER") sum(distances > config$caliper * sd(lp) + 1e-8) else NA_integer_
      if (nm == "CALIPER" && violations > 0) stop("卡尺匹配距离违反预设范围")
      retained <- weights[[nm]] > 0
      subclasses <- m$subclass[rownames(d)]
      cells <- table(subclasses[retained], factor(d$z[retained], levels = 0:1))
      if (any(rowSums(cells > 0) != 2)) stop("匹配子类缺少一个暴露组")
      matching_parts[[length(matching_parts) + 1L]] <- data.frame(imputation = i, method = nm,
        treated_retained = sum(retained & d$z == 1), control_retained = sum(retained & d$z == 0),
        subclasses = nrow(cells), caliper_violations = violations,
        max_pair_distance = if (length(distances)) max(distances) else NA_real_)
    }
  }
  for (nm in names(weights)) {
    w <- weights[[nm]]
    if (sum(w[d$z == 1]) == 0 || sum(w[d$z == 0]) == 0) stop("设计没有两组有效样本：", nm)
    balance_parts[[length(balance_parts) + 1L]] <- balance_one(d, p, w, design_estimands[[nm]]) |>
      mutate(imputation = i, design = nm, estimand = design_estimands[[nm]])
    weight_parts[[length(weight_parts) + 1L]] <- weight_stats(d, w) |>
      mutate(imputation = i, design = nm, estimand = design_estimands[[nm]])
  }
  design_objects[[i]] <- list(ps = pp, weights = weights, matches = matches,
    common_support = support, ps_decile_counts = decile_counts, clinical_cell_counts = clinical_counts)
}
positivity <- bind_rows(positivity_parts)
balance <- bind_rows(balance_parts)
weight_diagnostics <- bind_rows(weight_parts)
matching_quality <- bind_rows(matching_parts)
design_summary <- balance |>
  group_by(design, estimand, imputation) |>
  summarise(max_smd = max(abs_smd), max_ecdf = max(ecdf_difference), .groups = "drop") |>
  group_by(design, estimand) |>
  summarise(imputations = n(), median_max_smd = median(max_smd), worst_smd = max(max_smd),
    worst_ecdf = max(max_ecdf), all_balance_pass = all(max_smd < config$smd_limit), .groups = "drop") |>
  left_join(weight_diagnostics |> filter(group == "Overall") |>
    group_by(design) |> summarise(median_ess = median(ess), min_ess = min(ess),
      max_weight = max(max_weight), .groups = "drop"), by = "design")
matching_status <- design_summary |> filter(design %in% c("FULL", "NEAREST", "CALIPER", "GENETIC"))
# 仅所有插补都通过冻结平衡标准的匹配方案才可在下一阶段报告验证性效应。
# 失败方案保留诊断。算法增加不保证平衡，也不能用结局显著性选设计。
if (!all(design_summary$all_balance_pass[design_summary$design == "STABILIZED_IPTW"])) {
  stop("主要 IPTW 没有在所有插补中通过平衡门槛，请先重新评估设计")
}

if (config$draw_plots) {
  # 保留原项目的第 13 次代表插补；小规模验证用第 1 次，绝不挑最好看的一次。
  representative <- if (length(completed) >= config$representative) config$representative else 1L
  d <- completed[[representative]]; d$ps <- design_objects[[representative]]$ps$ps
  plots$ps <- ggplot2::ggplot(d, ggplot2::aes(ps, fill = qsmk)) +
    ggplot2::geom_density(alpha = .35) + ggplot2::labs(x = "Propensity score", fill = "Exposure") +
    ggplot2::theme_minimal()
  love_data <- balance |> filter(!feature %in% c("propensity_score", "logit_ps")) |>
    group_by(design, feature) |> summarise(abs_smd = max(abs_smd), .groups = "drop")
  plots$balance <- ggplot2::ggplot(love_data, ggplot2::aes(abs_smd, feature, colour = design)) +
    ggplot2::geom_point() + ggplot2::geom_vline(xintercept = config$smd_limit, linetype = 2) +
    ggplot2::labs(x = "Worst absolute SMD across imputations", y = NULL) + ggplot2::theme_minimal()
  if (config$display_plots) { print(plots$ps); print(plots$balance) }
}
print(design_summary)
