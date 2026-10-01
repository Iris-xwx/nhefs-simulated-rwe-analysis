# 1-5：缺失模式/选择模型、完整案例、MICE、Table 1/SMD/ECDF/共线性。
# 单独重跑本阶段前先运行 00_main.R 的参数和 functions.R 部分。
message("[1/3] 缺失与基线分析")
clean_raw <- readr::read_csv(file.path(data_dir, "nhefs_rwe_clean.csv"),
                           show_col_types = FALSE, col_types = readr::cols(id = readr::col_character()))
candidate_raw <- readr::read_csv(file.path(data_dir, "nhefs_rwe_analysis_candidate.csv"),
                               show_col_types = FALSE, col_types = readr::cols(id = readr::col_character()))
stopifnot(!anyDuplicated(candidate_raw$id), !anyDuplicated(clean_raw$id))
candidate <- typed_candidate(candidate_raw)
stopifnot(!anyNA(candidate$qsmk), !anyNA(candidate$outcome))
clean_input <- data.frame(id = clean_raw$id, source_system = clean_raw$source_system,
                         qsmk = clean_raw$qsmk_clean, outcome = clean_raw$wt82_71_clean)
for (v in baseline) clean_input[[v]] <- clean_raw[[paste0(v, "_clean")]]
clean <- typed_candidate(clean_input)
miss <- is.na(clean[, c("qsmk", "outcome", baseline)])
missingness <- data.frame(variable = colnames(miss), n_missing = colSums(miss),
                         percent = 100 * colMeans(miss), n_total = nrow(clean))
# 显式用列名标记组合，避免 c_across() 丢掉名称后把所有模式写成 Complete。
pattern <- apply(miss, 1, function(k) {
  if (any(k)) paste(colnames(miss)[k], collapse = ";") else "Complete"
})
missing_patterns <- as.data.frame(table(pattern), stringsAsFactors = FALSE)
names(missing_patterns) <- c("pattern", "n")
missing_patterns <- missing_patterns[order(-missing_patterns$n), ]
stratified_missing <- clean |>
  mutate(age_group = cut(age, c(-Inf, 39, 49, 59, Inf)),
    smoke_group = cut(smokeintensity, c(-Inf, 15, 25, Inf))) |>
  group_by(qsmk, source_system, age_group, smoke_group) |>
  summarise(n = n(), across(all_of(baseline), ~mean(is.na(.x))), .groups = "drop")

clean$any_missing <- rowSums(is.na(clean[baseline])) > 0
missing_model <- glm(any_missing ~ age + sex + race + source_system,
                     data = clean, family = binomial())
logistic_table <- function(fit) {
  cf <- coef(summary(fit)); ci <- confint.default(fit)
  data.frame(term = rownames(cf), odds_ratio = exp(cf[, 1]),
    conf_low = exp(ci[, 1]), conf_high = exp(ci[, 2]), p_value = cf[, 4])
}
missing_selection <- logistic_table(missing_model)
# 缺失与已观测变量有关支持 MCAR 不可信；不能据此证明 MAR，MI 仍依赖 MAR 工作假设。
candidate$cc <- complete.cases(candidate[baseline])
cc_data <- candidate[candidate$cc, ]
cc_model <- glm(cc ~ qsmk + source_system + age + sex + race,
                data = candidate, family = binomial())
cc_selection <- logistic_table(cc_model)
cc_population <- candidate |>
  group_by(cc) |> summarise(n = n(), quit_fraction = mean(qsmk == "Quit"),
    age_mean = mean(age, na.rm = TRUE), wt71_mean = mean(wt71, na.rm = TRUE), .groups = "drop")
cc_comparison <- tableone::CreateTableOne(vars = c("qsmk", baseline), strata = "cc",
  data = candidate, factorVars = c("qsmk", categorical), includeNA = TRUE, test = FALSE)
cc_comparison_output <- print(cc_comparison, smd = TRUE, printToggle = FALSE, showAllLevels = TRUE)

sample_flow <- data.frame(
  stage = c("Clean unique people", "Exposure known", "Exposure and outcome observed", "Analysis candidates", "CC", "MI"),
  n = c(nrow(clean), sum(!is.na(clean$qsmk)), sum(!is.na(clean$qsmk) & !is.na(clean$outcome)),
        nrow(candidate), nrow(cc_data), nrow(candidate))
)
tab1 <- tableone::CreateTableOne(vars = baseline, strata = "qsmk", data = cc_data,
                              factorVars = categorical, test = FALSE)
table1_output <- print(tab1, smd = TRUE, nonnormal = c("smokeintensity", "smokeyrs"),
                      showAllLevels = TRUE, printToggle = FALSE)
table1_smd <- tableone::ExtractSmd(tab1) # 此 Table 1 的人群是 CC；MI 后每次再检查平衡。

# 连续变量相关性、设计矩阵秩、标准化条件数及 VIF/GVIF。
# 临床重要的年龄/吸烟年限不因单个 VIF 较高而机械删除。
vif_fit <- lm(reformulate(baseline, response = "outcome"), data = cc_data)
x <- model.matrix(vif_fit)
collinearity <- list(
  correlation = cor(cc_data[continuous]), rank = qr(x)$rank, columns = ncol(x),
  scaled_condition_number = kappa(scale(x[, -1, drop = FALSE]), exact = TRUE),
  vif = if (qr(x)$rank == ncol(x)) car::vif(vif_fit) else alias(vif_fit)
)
# 模拟层中的 age_at_1982 与 wt71_lb_source 为确定性冗余，候选层已排除。

# MICE：只插补基线协变量；ID 仅作行名；观测的暴露/结局是预测变量。
if (anyNA(candidate$age)) stop("当前年龄平方辅助项要求 age 完整；年龄缺失时应改用被动插补")
mice_input <- candidate[, c("qsmk", "outcome", "source_system", "sex", "race", "age")]
mice_input$age_centered_sq <- (candidate$age - mean(candidate$age))^2
mice_input <- cbind(mice_input, candidate[, c("education", "smokeintensity", "smokeyrs", "exercise", "active", "wt71")])
initial <- mice::mice(mice_input, maxit = 0, printFlag = FALSE)
method <- setNames(rep("", ncol(mice_input)), names(mice_input))
to_impute <- baseline[vapply(candidate[baseline], anyNA, logical(1))]
method[intersect(to_impute, continuous)] <- "pmm"
method[intersect(to_impute, c("sex", "race"))] <- "logreg"
method[intersect(to_impute, c("education", "exercise", "active"))] <- "polyreg"
predictor <- matrix(1L, ncol(mice_input), ncol(mice_input),
                    dimnames = list(names(mice_input), names(mice_input)))
diag(predictor) <- 0L
predictor[method == "", ] <- 0L # 不插补变量仍可作为其他变量的预测列。
mice_methods <- data.frame(variable = names(method), method = unname(method),
                          n_missing = colSums(is.na(mice_input)))
imp <- mice::mice(mice_input, m = config$m, maxit = config$maxit,
                 method = method, predictorMatrix = predictor,
                 seed = config$seed, printFlag = FALSE)
completed <- lapply(seq_len(imp$m), function(i) analysis_data(mice::complete(imp, i)))
mice_events <- bind_rows(initial$loggedEvents, imp$loggedEvents)
mice_convergence <- mice::convergence(imp)
mice_postcheck <- bind_rows(lapply(seq_along(completed), function(i) {
  d <- completed[[i]]
  data.frame(imputation = i, missing_baseline = sum(is.na(d[baseline])),
    invalid_smoke = sum(d$smokeintensity < 0 | d$smokeintensity > 100),
    invalid_smokeyrs = sum(d$smokeyrs < 0 | d$smokeyrs > d$age),
    weight_review = sum(d$wt71 < 30 | d$wt71 > 180),
    unchanged_outcome = identical(d$outcome, candidate$outcome),
    unchanged_exposure = identical(d$qsmk, candidate$qsmk))
}))
stopifnot(all(mice_postcheck$missing_baseline == 0), all(mice_postcheck$unchanged_outcome),
          all(mice_postcheck$unchanged_exposure))
if (any(mice_postcheck[, c("invalid_smoke", "invalid_smokeyrs", "weight_review")] > 0)) {
  stop("插补后的范围/跨字段逻辑检查未通过，请检查插补模型，不自动更改插补值")
}
# 仅用预指定诊断模型量化有限 m 的模拟误差；其系数不是正式因果效应。
diagnostic <- with(imp, lm(outcome ~ qsmk + sex + race + age + education +
                          smokeintensity + smokeyrs + exercise + active + wt71))
diagnostic_pool <- mice::pool(diagnostic)$pooled
mice_mce <- diagnostic_pool |>
  mutate(mce = sqrt(b / m), mce_over_se = mce / sqrt(t)) |>
  select(term, m, mce, mce_over_se)
if (any(mice_mce$mce_over_se > .10)) warning("部分诊断系数 MCE/SE >0.10，建议增加 m")

plots <- list()
if (config$draw_plots) {
  plots$missing <- ggplot2::ggplot(missingness, ggplot2::aes(reorder(variable, percent), percent)) +
    ggplot2::geom_col(fill = "#387895") + ggplot2::coord_flip() +
    ggplot2::labs(x = NULL, y = "Missing (%)", title = "Missingness in the cleaned cohort") +
    ggplot2::theme_minimal()
  ecdf_long <- bind_rows(lapply(continuous, function(v) {
    data.frame(variable = v, value = cc_data[[v]], exposure = cc_data$qsmk)
  }))
  plots$ecdf <- ggplot2::ggplot(ecdf_long, ggplot2::aes(value, colour = exposure)) +
    ggplot2::stat_ecdf() + ggplot2::facet_wrap(~variable, scales = "free_x") +
    ggplot2::labs(y = "ECDF", title = "CC baseline distributions") + ggplot2::theme_minimal()
  plots$mice_trace <- plot(imp, y = to_impute, layout = c(2, length(to_impute))) # 全部变量同页。
  imputed_continuous <- intersect(to_impute, continuous)
  if (length(imputed_continuous)) plots$mice_density <- mice::densityplot(imp, reformulate(imputed_continuous))
  if (config$display_plots) for (p in plots) print(p)
}
