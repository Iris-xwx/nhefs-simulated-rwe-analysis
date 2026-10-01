# 公开教学数据 -> 人工模拟的 RWE 风格提取层 -> 确定性清洗。
# 用途：演示异构编码、缺失、单位与重复记录；不是新增的真实 EHR 研究。
# 从 00_main.R 调用。causaldata 提供源数据，dplyr 整理，readr 保存 CSV。
# 固定种子和随机调用顺序沿用外层原 00_prepare.R，不向原项目数据目录写入。
message("[0/3] 模拟质量降级与精简清洗")
data("nhefs", package = "causaldata", envir = environment())
source_nhefs <- as.data.frame(nhefs)
source_vars <- c("seqn", "qsmk", "wt82_71", baseline)
stopifnot(all(source_vars %in% names(source_nhefs)))

make_dirty_extract <- function(x, seed = 20260926L) {
  set.seed(seed)
  d <- x |> transmute(
    id = as.character(seqn), qsmk = if_else(qsmk == 1, "Quit", "Continued"),
    wt82_71 = as.numeric(wt82_71), sex = as.character(sex), race = as.character(race),
    age = as.numeric(age), education = as.numeric(education),
    smokeintensity = as.numeric(smokeintensity), smokeyrs = as.numeric(smokeyrs),
    exercise = as.character(exercise), active = as.character(active), wt71 = as.numeric(wt71)
  ) |> mutate(
    source_system = if_else(row_number() %% 4L == 0L, "Registry_B", "EHR_A"),
    extract_timestamp = if_else(source_system == "Registry_B", "2026-08-31 23:59:00", "2026-09-01 00:15:00"),
    age_at_1982 = age + 11, wt71_lb_source = wt71 * 2.2046226218
  )
  n <- nrow(d)
  # 1. 同义词、大小写、空格、前导零和无法解释的暴露编码。
  quit_rows <- which(d$qsmk == "Quit"); continue_rows <- which(d$qsmk == "Continued")
  d$qsmk[sample(quit_rows, 18L)] <- " quit "
  d$qsmk[sample(quit_rows, 9L)] <- "Y"
  d$qsmk[sample(continue_rows, 20L)] <- "CONTINUE"
  d$qsmk[sample(continue_rows, 8L)] <- "N"
  d$qsmk[sample(seq_len(n), 4L)] <- "Unknown"
  sex_rows <- sample(seq_len(n), 20L)
  d$sex[sex_rows[1:10]] <- paste0(" ", d$sex[sex_rows[1:10]])
  d$sex[sex_rows[11:20]] <- paste0(d$sex[sex_rows[11:20]], " ")
  race_rows <- sample(seq_len(n), 24L)
  d$race[race_rows[1:8]] <- paste0(" ", d$race[race_rows[1:8]])
  d$race[race_rows[9:16]] <- paste0(d$race[race_rows[9:16]], " ")
  d$race[race_rows[17:24]] <- paste0("0", d$race[race_rows[17:24]])
  # 2. MAR 风格缺失：依赖年龄、吸烟强度、人工来源标签；不证明真实缺失为 MAR。
  p_edu <- plogis(-3.3 + .035 * (d$age - mean(d$age, na.rm = TRUE)) +
                   .7 * (d$source_system == "Registry_B"))
  idx <- which(runif(n) < p_edu & !is.na(d$education)); d$education[idx] <- NA_real_
  p_exercise <- plogis(-3.2 + .018 * d$smokeintensity + .6 * (d$source_system == "Registry_B"))
  idx <- which(runif(n) < p_exercise & !is.na(d$exercise)); d$exercise[idx] <- NA_character_
  p_weight <- plogis(-3.6 + .5 * (d$source_system == "Registry_B") + .012 * (d$age - 45))
  idx <- which(runif(n) < p_weight & !is.na(d$wt71)); d$wt71[idx] <- NA_real_
  # 3. 哨兵、单位混用与不可能的吸烟年限；不改变已观测的结局。
  idx <- sample(which(!is.na(d$smokeintensity)), 7L); d$smokeintensity[idx] <- 999
  idx <- sample(which(!is.na(d$wt71)), 6L); d$wt71[idx] <- d$wt71[idx] * 2.2046226218
  idx <- sample(which(!is.na(d$age) & !is.na(d$smokeyrs)), 6L)
  d$smokeyrs[idx] <- d$age[idx] + sample(1:4, 6L, replace = TRUE)
  # 4. 重复提取：晚到版本不应机械覆盖更完整的基线记录。
  dup_idx <- sample(seq_len(n), 12L); d_dup <- d[dup_idx, , drop = FALSE]
  d_dup$extract_timestamp <- "2026-09-02 08:30:00"; d_dup$source_system <- "Late_feed"
  d_dup$exercise[1:3] <- NA_character_
  bind_rows(d, d_dup) |> mutate(extract_row_id = sprintf("ROW%05d", row_number())) |>
    select(extract_row_id, everything())
}

dirty_path <- file.path(data_dir, "nhefs_rwe_dirty.csv")
if (config$rebuild || !file.exists(dirty_path)) {
  readr::write_csv(make_dirty_extract(source_nhefs, config$preparation_seed), dirty_path, na = "")
}
# 原始 CSV 按字符读，保留空格与前导零；清洗不读取源数据真值恢复被删除的基线。
dirty <- readr::read_csv(dirty_path, col_types = readr::cols(.default = readr::col_character()),
                       na = c("", "NA"), trim_ws = FALSE, show_col_types = FALSE)
stopifnot(!anyNA(dirty$id), !anyDuplicated(dirty$extract_row_id))
number <- function(x) suppressWarnings(as.numeric(trimws(x)))
category <- function(x, allowed) { x <- number(x); ifelse(x %in% allowed, x, NA_real_) }
lb_per_kg <- 2.2046226218
stage <- dirty
stage$qsmk_clean <- case_when(tolower(trimws(stage$qsmk)) %in% c("quit", "y") ~ 1L,
  tolower(trimws(stage$qsmk)) %in% c("continued", "continue", "n") ~ 0L, TRUE ~ NA_integer_)
for (v in c("sex", "race")) stage[[paste0(v, "_clean")]] <- category(stage[[v]], 0:1)
stage$education_clean <- category(stage$education, 1:5)
for (v in c("exercise", "active")) stage[[paste0(v, "_clean")]] <- category(stage[[v]], 0:2)
age <- number(stage$age); smoke <- number(stage$smokeintensity); years <- number(stage$smokeyrs)
stage$age_clean <- ifelse(age >= 18 & age <= 100, age, NA_real_)
stage$smokeintensity_clean <- ifelse(smoke >= 0 & smoke <= 100 & smoke != 999, smoke, NA_real_)
stage$smokeyrs_clean <- ifelse(years >= 0 & years <= age, years, NA_real_)
w <- number(stage$wt71); reference <- number(stage$wt71_lb_source)
kg_match <- !is.na(w) & !is.na(reference) & abs(w - reference/lb_per_kg) <= 1
lb_match <- !is.na(w) & !is.na(reference) & abs(w - reference) <= 1
unresolved <- !is.na(w) & !is.na(reference) & !kg_match & !lb_match
stage$wt71_clean <- ifelse(unresolved, NA_real_, ifelse(lb_match, w/lb_per_kg, w))
stage$weight_unit_status <- case_when(is.na(w) ~ "SOURCE_MISSING", lb_match ~ "LB_VALUE_CONVERTED_TO_KG",
  kg_match ~ "KG_CONFIRMED", unresolved ~ "UNIT_UNRESOLVED", TRUE ~ "REFERENCE_MISSING")
# wt71_lb_source 是模拟时产生的单位谱系，并非独立真实校验；不凭体重大小猜单位。
stage$wt82_71_clean <- number(stage$wt82_71)
stage$review_required <- unresolved | (!is.na(stage$wt71_clean) &
                                       (stage$wt71_clean < 30 | stage$wt71_clean > 180))
fields <- paste0(baseline, "_clean")
stage$baseline_complete_count <- rowSums(!is.na(stage[fields]))
stage$source_priority <- unname(c(EHR_A = 1L, Registry_B = 2L, Late_feed = 3L)[stage$source_system])
if (anyNA(stage$source_priority)) stop("存在未定义来源：请预先定义版本优先级")
stage$timestamp <- as.POSIXct(stage$extract_timestamp, format = "%Y-%m-%d %H:%M:%S", tz = "UTC")
stopifnot(!anyNA(stage$timestamp))
ranked <- stage |> group_by(id) |>
  arrange(desc(baseline_complete_count), source_priority, desc(timestamp), extract_row_id, .by_group = TRUE) |>
  mutate(version_rank = row_number()) |> ungroup()
clean_export <- ranked |> filter(version_rank == 1L) |> select(-timestamp, -version_rank) |>
  mutate(any_baseline_missing = if_any(all_of(fields), is.na), analysis_status = case_when(
    is.na(qsmk_clean) ~ "EXCLUDED_EXPOSURE_UNRESOLVED",
    is.na(wt82_71_clean) ~ "EXCLUDED_OUTCOME_MISSING", review_required ~ "REVIEW_REQUIRED",
    any_baseline_missing ~ "ELIGIBLE_FOR_BASELINE_MI", TRUE ~ "COMPLETE_CASE_ELIGIBLE"))
quarantine <- ranked |> filter(version_rank > 1L) |> select(-timestamp) |>
  mutate(quarantine_reason = "NON_SELECTED_DUPLICATE_VERSION")
candidate_export <- clean_export |> filter(analysis_status %in%
  c("ELIGIBLE_FOR_BASELINE_MI", "COMPLETE_CASE_ELIGIBLE")) |>
  transmute(id, source_system, extract_timestamp, qsmk = qsmk_clean, outcome = wt82_71_clean,
    sex = sex_clean, race = race_clean, age = age_clean, education = education_clean,
    smokeintensity = smokeintensity_clean, smokeyrs = smokeyrs_clean,
    exercise = exercise_clean, active = active_clean, wt71 = wt71_clean)
stopifnot(nrow(clean_export) + nrow(quarantine) == nrow(dirty), !anyDuplicated(clean_export$id),
          !anyNA(candidate_export$qsmk), !anyNA(candidate_export$outcome))
for (nm in c("clean", "analysis_candidate", "quarantine")) {
  obj <- switch(nm, clean = clean_export, analysis_candidate = candidate_export, quarantine = quarantine)
  readr::write_csv(obj, file.path(data_dir, paste0("nhefs_rwe_", nm, ".csv")), na = "")
}
# 仅为报告比较缺失来源；真值不进入清洗或插补模型。
base_rows <- dirty[match(as.character(source_nhefs$seqn), dirty$id), ]
simulation_summary <- bind_rows(lapply(c("qsmk", "wt82_71", baseline), function(v) {
  raw <- base_rows[[v]]
  changed_na <- if (v == "qsmk") raw == "Unknown" else is.na(raw)
  data.frame(variable = v, original_missing = sum(is.na(source_nhefs[[v]])),
    dirty_missing_or_unknown = sum(changed_na, na.rm = TRUE),
    newly_missing_or_unknown = sum(changed_na & !is.na(source_nhefs[[v]]), na.rm = TRUE))
}))
cleaning_summary <- data.frame(
  issue = c("非标准但可映射的戒烟编码", "不可解析的暴露", "每日吸烟量哨兵 999", "基线体重磅值",
            "吸烟年限大于年龄", "重复提取行", "清洗后待复核", "结局原有缺失"),
  affected_rows = c(sum(!is.na(stage$qsmk_clean) & !stage$qsmk %in% c("Quit", "Continued")),
    sum(is.na(stage$qsmk_clean)), sum(smoke == 999, na.rm = TRUE), sum(lb_match),
    sum(years > age, na.rm = TRUE), nrow(quarantine), sum(clean_export$review_required),
    sum(is.na(source_nhefs$wt82_71))),
  counting_population = c(rep("dirty 提取行；去重前", 5), "隔离提取行", "去重后患者", "原公开数据患者"),
  action = c("词典标准化", "排除主分析；不推测暴露", "置 NA，交给基线 MICE",
    "按已知谱系换算 kg", "置 NA，不截成合法值", "完整性→来源→时间→行号排序保留一个版本",
    "未经复核不进入主分析", "主分析排除，不擅自插补结局")
)
preparation_flow <- data.frame(stage = c("Public source people", "Dirty extract rows", "Quarantined duplicate rows"),
                               n = c(nrow(source_nhefs), nrow(dirty), nrow(quarantine)))
