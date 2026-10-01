# 手动执行的安装助手；不会在 00_main.R 中静默安装。
# 安装 CRAN 当前版本；严格版本复现另用 renv.lock（见 README）。
# causaldata 源数据；dplyr/readr 数据操作；ggplot2 科学图；mice 插补/Rubin；
# tableone 基线描述；car VIF；sandwich 稳健方差；MatchIt 四种匹配；
# optmatch/Matching/rgenoud 匹配后端；knitr/rmarkdown 供 Quarto 执行 R 报告。
packages <- c("causaldata", "dplyr", "readr", "ggplot2", "mice", "tableone", "car",
  "sandwich", "MatchIt", "optmatch", "Matching", "rgenoud", "knitr", "rmarkdown")
missing <- packages[!vapply(packages, requireNamespace, logical(1), quietly = TRUE)]
if (length(missing)) install.packages(missing, repos = "https://cloud.r-project.org")
message("依赖检查完成。splines 为 R 推荐包；Quarto 应用需单独安装。")
