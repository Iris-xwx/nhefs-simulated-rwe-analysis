# 上传前的只读检查。readr 读取结果；不重跑模型、不联网、不初始化 Git。
root <- if (file.exists("00_main.R")) "." else "nhefs"
root <- normalizePath(root, winslash = "/", mustWork = TRUE)
source(file.path(root, "functions.R"), encoding = "UTF-8")
read_result <- function(nm) readr::read_csv(file.path(root, "results", paste0(nm, ".csv")), show_col_types = FALSE)
hash <- read_result("code_checksums")
stopifnot(identical(unname(code_md5(file.path(root, hash$file))), hash$md5))
input_hash <- read_result("input_checksums")
available <- file.exists(file.path(root, "data", input_hash$file))
if (any(available)) stopifnot(identical(unname(code_md5(file.path(root, "data", input_hash$file[available]))), input_hash$md5[available]))
meta <- read_result("run_metadata"); value <- setNames(meta$value, meta$setting)
stopifnot(value[["m"]] == "20", value[["maxit"]] == "20", value[["bootstrap"]] == "500",
          value[["run_matching"]] == "TRUE", value[["draw_plots"]] == "TRUE")
flow <- read_result("sample_flow"); sizes <- setNames(flow$n, flow$stage)
stopifnot(sizes[["Dirty extract rows"]] == sizes[["Clean unique people"]] + sizes[["Quarantined duplicate rows"]],
          sizes[["CC"]] <= sizes[["MI"]], sizes[["MI"]] == sizes[["Analysis candidates"]])
effects <- read_result("pooled_effects")
stopifnot(all(is.finite(effects$estimate)), all(effects$se > 0),
          all(effects$conf_low <= effects$estimate & effects$conf_high >= effects$estimate))
primary <- effects[effects$method %in% c("IPTW", "AIPW") & effects$estimand == "ATE", ]
stopifnot(nrow(primary) == 2L, all(primary$m == 20L))
mi <- read_result("mice_diagnostics")
stopifnot(nrow(mi) == 20L, all(mi$missing_baseline == 0), all(mi$unchanged_outcome), all(mi$unchanged_exposure))
quality <- read_result("matching_quality")
stopifnot(nrow(quality) == 80L, all(quality$caliper_violations[quality$method == "CALIPER"] == 0L))
design <- read_result("design_summary")
stopifnot(design$all_balance_pass[design$design == "STABILIZED_IPTW"])
boot <- read_result("bootstrap")
stopifnot(nrow(boot) == 2L, all(boot$requested == 500), all(boot$success >= 475))
figures <- c("missing", "ecdf", "mice_trace", "mice_density", "ps", "balance", "effects")
stopifnot(all(file.exists(file.path(root, "results", "figures", paste0(figures, ".png")))),
          file.exists(file.path(root, "report.html")), file.exists(file.path(root, "renv.lock")))
# 可选 xml2：检查已渲染的正文与内嵌图，而不是误把隐藏的报告源码当作正文。
if (requireNamespace("xml2", quietly = TRUE)) {
  doc <- xml2::read_html(file.path(root, "report.html"))
  main <- xml2::xml_find_first(doc, "//main")
  body_text <- xml2::xml_text(main)
  image_paths <- xml2::xml_attr(xml2::xml_find_all(main, ".//img"), "src")
  stopifnot(grepl("正式默认参数运行", body_text, fixed = TRUE),
    !grepl("<U+", body_text, fixed = TRUE),
    length(image_paths) == 7L, all(startsWith(image_paths, "data:image/png;base64,")))
}
public_csv <- list.files(file.path(root, "results"), pattern = "[.]csv$", full.names = TRUE)
for (f in public_csv) {
  columns <- names(readr::read_csv(f, show_col_types = FALSE))
  stopifnot(!any(columns %in% c("id", "seqn", "most_influential_id")))
}
message("检查通过：正式参数、源码/输入、样本守恒、MICE、匹配、效应、bootstrap、图及报告。")
