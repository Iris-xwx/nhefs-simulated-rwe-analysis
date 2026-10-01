# 4、12-16：Rubin、标准化、AIPW、sandwich/IF 方差、bootstrap 和核心敏感性。
# 主分析：候选总体 ATE / 稳定 IPTW；AIPW 是核心双重稳健参照。
# ATT、Overlap、PS 修剪后的总体、权重截尾后的目标分别标记，不能当成同一个 ATE。
message("[3/3] 效应估计与敏感性分析")

estimate_core <- function(d, ps = NULL) {
  if (is.null(ps)) ps <- fit_ps(d, formulas$ps2)$ps
  w <- make_weights(d$z, ps, config$truncation)
  o <- fit_outcome(d)
  a <- aipw_effect(d, ps, o)
  bind_rows(
    weighted_effect(d, w$STABILIZED_IPTW, "IPTW"), a$row,
    or_standardized(d, o), or_standardized(d, flexible = FALSE),
    weighted_effect(d, w$UNADJUSTED, "UNADJUSTED", "Descriptive"),
    weighted_effect(d, w$TRUNCATED_IPTW, "TRUNCATED_IPTW", "Modified-weight ATE"),
    weighted_effect(d, w$ATT, "ATT_WEIGHTING", "ATT"),
    weighted_effect(d, w$OVERLAP, "OVERLAP", "Overlap")
  )
}
effect_parts <- influence_parts <- population_parts <- list()
trim_membership <- list()
passing_matches <- matching_status$design[matching_status$all_balance_pass]
for (i in seq_along(completed)) {
  d <- completed[[i]]; obj <- design_objects[[i]]
  p <- obj$ps$ps
  effect_parts[[length(effect_parts) + 1L]] <- estimate_core(d, p) |> mutate(imputation = i)
  a <- aipw_effect(d, p, fit_outcome(d))
  influence_parts[[i]] <- a$influence |> mutate(imputation = i)
  # 模型规格敏感性：线性 PS-1，OR-2 保持固定。
  p1 <- fit_ps(d, formulas$ps1)$ps
  w1 <- make_weights(d$z, p1, config$truncation)$STABILIZED_IPTW
  effect_parts[[length(effect_parts) + 1L]] <- bind_rows(
    weighted_effect(d, w1, "IPTW_PS1"), aipw_effect(d, p1, fit_outcome(d), "AIPW_PS1")$row
  ) |> mutate(imputation = i)

  # PS 修剪删除患者：0.05/0.95 和经验共同支持域分别保留，随后重新拟合 PS。
  # 阈值和模型提前固定；各插补成员变化意味着目标人群也可能略有变化。
  keep_sets <- list(
    PS_TRIM = p >= config$ps_trim[1] & p <= config$ps_trim[2],
    COMMON_SUPPORT = p >= obj$common_support[1] & p <= obj$common_support[2]
  )
  for (nm in names(keep_sets)) {
    keep <- keep_sets[[nm]]
    sub <- d[keep, ]
    sub_effects <- estimate_core(sub) |> filter(method %in% c("IPTW", "AIPW")) |>
      mutate(method = paste(nm, method, sep = "_"), estimand = paste(nm, "restricted ATE"), imputation = i)
    effect_parts[[length(effect_parts) + 1L]] <- sub_effects
    trim_membership[[length(trim_membership) + 1L]] <- data.frame(imputation = i,
      scenario = nm, id = d$id, included = keep)
    population_parts[[length(population_parts) + 1L]] <- data.frame(imputation = i,
      scenario = nm, n = nrow(sub), quit_fraction = mean(sub$z),
      age_mean = mean(sub$age), smoke_mean = mean(sub$smokeintensity), wt71_mean = mean(sub$wt71))
  }
  # ATT 和 overlap 的目标协变量分布也要展示，避免数值相近被解释成验证 ATE。
  for (nm in c("ATT", "OVERLAP")) {
    w <- obj$weights[[nm]]
    population_parts[[length(population_parts) + 1L]] <- data.frame(imputation = i,
      scenario = nm, n = nrow(d), quit_fraction = weighted.mean(d$z, w),
      age_mean = weighted.mean(d$age, w), smoke_mean = weighted.mean(d$smokeintensity, w),
      wt71_mean = weighted.mean(d$wt71, w))
  }
  # 匹配的正式效应仅在设计通过门槛后估计；variance 按 subclass 聚类。
  for (nm in passing_matches) {
    effect_parts[[length(effect_parts) + 1L]] <- weighted_effect(
      d, obj$weights[[nm]], paste0("MATCH_", nm), design_estimands[[nm]],
      subclass = as.character(obj$matches[[nm]]$subclass[rownames(d)])
    ) |> mutate(imputation = i)
  }
}

# 完整案例是不同缺失处理/选入人群敏感性：模型结点仍用预先冻结的规格。
cc_effect <- estimate_core(analysis_data(cc_data)) |>
  mutate(method = paste0("CC_", method), estimand = paste("CC", estimand), imputation = 0L)
effect_by_imputation <- bind_rows(effect_parts, list(cc_effect))
pooled_effects <- pool_effects(effect_by_imputation)
influence_values <- bind_rows(influence_parts)
influence_diagnostics <- influence_values |> group_by(imputation) |>
  summarise(max_abs = max(abs(influence)),
    most_influential_id = id[which.max(abs(influence))],
    max_absolute_share = max(abs(influence)) / sum(abs(influence)),
    top_1pct_share = sum(head(sort(abs(influence), decreasing = TRUE),
      ceiling(n() * .01))) / sum(abs(influence)), .groups = "drop")
sensitivity_populations <- bind_rows(population_parts)
trim_inclusion <- bind_rows(trim_membership) |> group_by(scenario, id) |>
  summarise(inclusion_fraction = mean(included), .groups = "drop")

# Bootstrap 抽样单位是参与者；每次从完成数据重抽样后重新拟合 PS-2 和 OR-2。
# 固定一个代表插补，结果只作为方差敏感性；不声称重复了 MICE 全流程。
bootstrap_results <- data.frame()
bootstrap_replicates <- data.frame()
if (config$bootstrap > 0L) {
  representative <- if (length(completed) >= config$representative) config$representative else 1L
  boot_data <- completed[[representative]]
  set.seed(config$bootstrap_seed)
  boot_list <- lapply(seq_len(config$bootstrap), function(b) {
    if (b %% 100L == 0L) message("  bootstrap ", b, "/", config$bootstrap)
    d <- boot_data[sample.int(nrow(boot_data), replace = TRUE), ]
    tryCatch({
      p <- fit_ps(d, formulas$ps2)$ps
      w <- make_weights(d$z, p, config$truncation)$STABILIZED_IPTW
      data.frame(replicate = b, iptw = weighted_effect(d, w, "IPTW")$estimate,
        aipw = aipw_effect(d, p, fit_outcome(d))$row$estimate, error = NA_character_)
    }, error = function(e) data.frame(replicate = b, iptw = NA_real_, aipw = NA_real_,
                                     error = conditionMessage(e)))
  })
  bootstrap_replicates <- bind_rows(boot_list)
  bootstrap_results <- bind_rows(lapply(c("iptw", "aipw"), function(v) {
    x <- bootstrap_replicates[[v]]; x <- x[is.finite(x)]
    if (length(x) < 2L) stop("有效 bootstrap 少于两次：", v)
    original <- effect_by_imputation |>
      filter(imputation == representative, method == toupper(v), estimand == "ATE")
    data.frame(method = toupper(v), representative_imputation = representative,
      success = length(x), requested = config$bootstrap, original_estimate = original$estimate,
      analytic_se = sqrt(original$variance), bootstrap_se = sd(x),
      se_ratio = sd(x) / sqrt(original$variance), percentile_low = unname(quantile(x, .025)),
      percentile_high = unname(quantile(x, .975)))
  }))
  if (any(bootstrap_results$success < .95 * config$bootstrap)) warning("bootstrap 失败率超过 5%，请检查失败原因")
}

if (config$draw_plots) {
  # 按目标人群分面展示；不把 ATT/overlap/CC/修剪估计混作同一总体的验证结果。
  effects_plot_data <- pooled_effects |> mutate(
    display_method = sub("^(CC_|COMMON_SUPPORT_|PS_TRIM_)", "", method),
    display_target = case_when(estimand == "COMMON_SUPPORT restricted ATE" ~ "Common-support ATE",
      estimand == "PS_TRIM restricted ATE" ~ "PS-trimmed ATE", TRUE ~ estimand))
  plots$effects <- ggplot2::ggplot(effects_plot_data,
    ggplot2::aes(estimate, reorder(display_method, estimate), xmin = conf_low, xmax = conf_high)) +
    ggplot2::geom_pointrange() + ggplot2::geom_vline(xintercept = 0, linetype = 2) +
    ggplot2::facet_wrap(~display_target, scales = "free_y", ncol = 2,
                        labeller = ggplot2::label_wrap_gen(width = 26)) +
    ggplot2::labs(x = "Mean difference (kg), 95% CI", y = NULL) + ggplot2::theme_minimal()
  if (config$display_plots) print(plots$effects)
}
# 解释仍依赖一致性、条件交换性和正值性。当前脚本覆盖图中主体方法；
# 不执行完整版本额外的 E-value、结局 MAR/MNAR、entropy calibration 与交叉拟合。
# 结局缺失、戒烟时点与 time zero 不完全对齐、治疗版本和外部可推广性仍为局限。
