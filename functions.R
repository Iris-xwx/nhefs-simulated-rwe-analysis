# 共用方法函数：避免在每个分析阶段重复 PS 公式、平衡计算和效应估计。
baseline <- c("sex", "race", "age", "education", "smokeintensity",
              "smokeyrs", "exercise", "active", "wt71")
continuous <- c("age", "smokeintensity", "smokeyrs", "wt71")
categorical <- setdiff(baseline, continuous)

# 报告一致性检查：统一 UTF-8/LF 后计算指纹，避免 Git 换行转换导致误报。
code_md5 <- function(paths) {
  vapply(paths, function(p) {
    tmp <- tempfile(fileext = ".txt")
    on.exit(unlink(tmp), add = TRUE)
    lines <- readLines(p, encoding = "UTF-8", warn = FALSE)
    writeBin(charToRaw(enc2utf8(paste0(paste(lines, collapse = "\n"), "\n"))), tmp)
    unname(tools::md5sum(tmp))
  }, character(1))
}

typed_candidate <- function(x) {
  data.frame(
    id = as.character(x$id), source_system = factor(x$source_system),
    qsmk = factor(x$qsmk, 0:1, c("Continued", "Quit")), outcome = x$outcome,
    sex = factor(x$sex, 0:1, c("Male", "Female")),
    race = factor(x$race, 0:1, c("White", "Black_or_other")), age = x$age,
    education = factor(x$education, 1:5, c("8th_grade_or_less", "High_school_dropout",
      "High_school", "College_dropout", "College_or_higher")),
    smokeintensity = x$smokeintensity, smokeyrs = x$smokeyrs,
    exercise = factor(x$exercise, 0:2, c("Much", "Moderate", "No")),
    active = factor(x$active, 0:2, c("Very_active", "Moderately_active", "Inactive")),
    wt71 = x$wt71, row.names = as.character(x$id)
  )
}
analysis_data <- function(d) {
  d$z <- as.integer(d$qsmk == "Quit")
  d$id <- rownames(d)
  d
}
ess <- function(w) {
  w <- w[is.finite(w) & w > 0]
  if (!length(w)) return(NA_real_)
  sum(w)^2 / sum(w^2)
}

# 模型结点来自未插补候选数据的观测值，跨插补、CC、修剪和 bootstrap 固定。
make_formulas <- function(candidate) {
  knots <- lapply(candidate[continuous], function(x) {
    x <- x[is.finite(x)]
    list(knots = unname(quantile(x, c(1/3, 2/3))), boundary = range(x))
  })
  ns_terms <- vapply(continuous, function(v) {
    s <- knots[[v]]
    sprintf("splines::ns(%s, knots=c(%s), Boundary.knots=c(%s))", v,
      paste(s$knots, collapse = ","), paste(s$boundary, collapse = ","))
  }, character(1))
  linear_rhs <- paste(baseline, collapse = " + ")
  flexible_rhs <- paste(c(categorical, ns_terms), collapse = " + ")
  list(
    ps1 = as.formula(paste("z ~", linear_rhs)),
    ps2 = as.formula(paste("z ~", flexible_rhs)),
    or1 = as.formula(paste("outcome ~ z +", linear_rhs)),
    or2 = as.formula(paste("outcome ~", flexible_rhs)), knots = knots
  )
}
fit_ps <- function(d, formula) {
  fit <- glm(formula, family = binomial(), data = d)
  if (!fit$converged || anyNA(coef(fit))) stop("PS Logistic 拟合失败或设计矩阵不满秩")
  raw_ps <- as.numeric(predict(fit, type = "response"))
  # clip 仅防止数值除零；未经 clip 的分布仍用于正值性诊断。
  list(fit = fit, raw_ps = raw_ps, ps = pmin(pmax(raw_ps, 1e-6), 1 - 1e-6))
}

make_weights <- function(z, p, truncation) {
  pr <- mean(z)
  sw <- ifelse(z == 1, pr / p, (1 - pr) / (1 - p))
  limits <- quantile(sw, truncation)
  list(
    UNADJUSTED = rep(1, length(z)),
    IPTW = ifelse(z == 1, 1 / p, 1 / (1 - p)),
    STABILIZED_IPTW = sw,
    ATT = ifelse(z == 1, 1, p / (1 - p)),
    OVERLAP = ifelse(z == 1, 1 - p, p),
    TRUNCATED_IPTW = pmin(pmax(sw, limits[1]), limits[2])
  )
}
support_limits <- function(z, p) c(
  max(min(p[z == 0]), min(p[z == 1])),
  min(max(p[z == 0]), max(p[z == 1]))
)

# 固定使用调整前的标准差作 SMD 分母，保证各设计之间可比较。
# ATT 分母用戒烟者 SD；ATE/overlap 用两组调整前合并 SD。
# ECDF 差比较整个分布；单个分类水平也作为 0/1 特征检查。
balance_features <- function(d, ps) {
  x <- as.matrix(d[continuous])
  for (v in categorical) {
    for (level in levels(d[[v]])) {
      value <- matrix(as.numeric(d[[v]] == level), ncol = 1,
                      dimnames = list(NULL, paste(v, level, sep = "::")))
      x <- cbind(x, value) # 包括参考类别，避免漏掉多分类变量最失衡的水平。
    }
  }
  for (v in continuous) {
    s <- formulas$knots[[v]]
    basis <- splines::ns(d[[v]], knots = s$knots, Boundary.knots = s$boundary)
    colnames(basis) <- paste0(v, "::ns", seq_len(ncol(basis)))
    x <- cbind(x, basis)
  }
  cbind(x, propensity_score = ps, logit_ps = qlogis(ps))
}
weighted_ecdf <- function(x, w, grid) {
  if (sum(w) <= 0) return(rep(NA_real_, length(grid)))
  o <- order(x); k <- findInterval(grid, x[o]); cw <- c(0, cumsum(w[o]) / sum(w))
  cw[k + 1L]
}
balance_one <- function(d, ps, w, target = "ATE") {
  x <- balance_features(d, ps)
  z <- d$z
  bind_rows(lapply(seq_len(ncol(x)), function(j) {
    v <- x[, j]; v1 <- v[z == 1]; v0 <- v[z == 0]
    w1 <- w[z == 1]; w0 <- w[z == 0]
    denominator <- if (target == "ATT") sd(v1) else sqrt((var(v1) + var(v0))/2)
    difference <- weighted.mean(v1, w1) - weighted.mean(v0, w0)
    smd <- if (denominator > 0) difference / denominator else if (abs(difference) < 1e-12) 0 else Inf
    grid <- sort(unique(v))
    data.frame(feature = colnames(x)[j], smd = smd, abs_smd = abs(smd),
      ecdf_difference = max(abs(weighted_ecdf(v1, w1, grid) - weighted_ecdf(v0, w0, grid))))
  }))
}
weight_stats <- function(d, w) {
  bind_rows(lapply(c("Overall", "Quit", "Continued"), function(g) {
    k <- if (g == "Overall") rep(TRUE, nrow(d)) else d$qsmk == g
    v <- w[k]
    data.frame(group = g, n = sum(k), retained = sum(v > 0),
      ess = ess(v), max_weight = max(v), p99_weight = unname(quantile(v, .99)))
  }))
}

# 多重插补：先在每个数据集中估计 Q 和 U，再做 Rubin 合并，不堆叠或平均 PS。
# pool.scalar 采用有限样本 Barnard-Rubin 自由度，并给出真正的 FMI。
pool_effects <- function(rows) {
  keys <- interaction(rows$method, rows$estimand, drop = TRUE)
  bind_rows(lapply(split(rows, keys), function(x) {
    if (nrow(x) == 1L) {
      q <- x$estimate; se <- sqrt(x$variance); df <- max(1, x$n - x$parameters)
      b <- 0; u <- x$variance; fmi <- NA_real_
    } else {
      p <- mice::pool.scalar(x$estimate, x$variance, n = min(x$n), k = max(x$parameters))
      q <- p$qbar; se <- sqrt(p$t); df <- p$df; b <- p$b; u <- p$ubar; fmi <- p$fmi
    }
    data.frame(method = x$method[1], estimand = x$estimand[1], m = nrow(x),
      estimate = q, se = se, conf_low = q - qt(.975, df)*se,
      conf_high = q + qt(.975, df)*se, df = df, between_variance = b,
      within_variance = u, fmi = fmi)
  }))
}
effect_row <- function(method, estimand, q, u, n, parameters = 1L) {
  if (!is.finite(q) || !is.finite(u) || u < 0) stop("效应估计或方差无效：", method)
  data.frame(method = method, estimand = estimand, estimate = q,
    variance = u, n = n, parameters = parameters)
}
weighted_effect <- function(d, w, method, target = "ATE", subclass = NULL) {
  keep <- w > 0
  if (any(table(factor(d$z[keep], levels = 0:1)) == 0)) stop("加权后缺少一个暴露组")
  x <- d[keep, ]; x$.w <- w[keep]
  fit <- lm(outcome ~ z, data = x, weights = .w)
  # 非匹配估计视权重为固定；匹配使用 subclass 聚类，不能当成独立观察。
  vc <- if (is.null(subclass)) sandwich::vcovHC(fit, type = "HC0") else {
    cl <- subclass[keep]
    if (anyNA(cl) || length(unique(cl)) < 2L) stop("匹配子类不足以估计聚类方差")
    sandwich::vcovCL(fit, cluster = cl, type = "HC0")
  }
  row <- effect_row(method, target, coef(fit)[["z"]], vc["z", "z"], nrow(x), 2L)
  row$mean_quit <- weighted.mean(x$outcome[x$z == 1], x$.w[x$z == 1])
  row$mean_continued <- weighted.mean(x$outcome[x$z == 0], x$.w[x$z == 0])
  row
}
fit_outcome <- function(d) {
  fit1 <- lm(formulas$or2, d[d$z == 1, ])
  fit0 <- lm(formulas$or2, d[d$z == 0, ])
  if (anyNA(coef(fit1)) || anyNA(coef(fit0))) stop("OR-2 模型不满秩")
  list(fit1 = fit1, fit0 = fit0,
    m1 = as.numeric(predict(fit1, d)), m0 = as.numeric(predict(fit0, d)))
}
or_standardized <- function(d, outcome_fit = NULL, flexible = TRUE) {
  if (!flexible) {
    fit <- lm(formulas$or1, d)
    d1 <- d0 <- d; d1$z <- 1; d0$z <- 0
    tt <- delete.response(terms(fit))
    x1 <- model.matrix(tt, d1); x0 <- model.matrix(tt, d0)
    g <- colMeans(x1 - x0)
    q <- sum(g * coef(fit)); u <- drop(t(g) %*% sandwich::vcovHC(fit, type = "HC0") %*% g)
    row <- effect_row("OR1", "ATE", q, u, nrow(d), length(coef(fit)))
    row$mean_quit <- mean(predict(fit, d1)); row$mean_continued <- mean(predict(fit, d0))
    return(row)
  }
  o <- outcome_fit
  g1 <- colMeans(model.matrix(delete.response(terms(o$fit1)), d))
  g0 <- colMeans(model.matrix(delete.response(terms(o$fit0)), d))
  u <- drop(t(g1) %*% sandwich::vcovHC(o$fit1, type = "HC0") %*% g1 +
              t(g0) %*% sandwich::vcovHC(o$fit0, type = "HC0") %*% g0)
  # delta 方差条件于当前经验协变量分布；随机目标人群推断以 AIPW/重拟合 bootstrap 为参照。
  row <- effect_row("OR2", "ATE", mean(o$m1 - o$m0), u, nrow(d),
    max(length(coef(o$fit1)), length(coef(o$fit0))))
  row$mean_quit <- mean(o$m1); row$mean_continued <- mean(o$m0)
  row
}

aipw_effect <- function(d, ps, o, method = "AIPW") {
  phi1 <- o$m1 + d$z * (d$outcome - o$m1) / ps
  phi0 <- o$m0 + (1 - d$z) * (d$outcome - o$m0) / (1 - ps)
  phi <- phi1 - phi0
  row <- effect_row(method, "ATE", mean(phi), var(phi) / nrow(d), nrow(d))
  row$mean_quit <- mean(phi1); row$mean_continued <- mean(phi0)
  # 同样本拟合的经验影响函数方差是近似；bootstrap 重拟合 PS 和 OR 检查其稳健性。
  list(row = row, influence = data.frame(id = d$id, influence = phi - mean(phi)))
}
