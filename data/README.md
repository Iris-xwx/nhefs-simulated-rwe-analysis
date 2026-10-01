# 数据来源与再生成

源数据是 R 包 `causaldata` 中的公开教学数据 `nhefs`，由 NHEFS 观察性随访数据整理而来；不是本项目采集的临床 EHR。

- [Hernán 与 Robins 的教材及数据](https://miguelhernan.org/whatifbook)
- [causaldata 包及源代码](https://github.com/NickCH-K/causaldata)

运行 `source("00_main.R", encoding="UTF-8")` 后，本目录生成模拟提取表、清洗表、分析候选表和重复版本隔离表。不需要外层项目的 `data/`。

`EHR_A`、`Registry_B`、`Late_feed`、提取时间、附加缺失、哨兵值、混合单位和重复记录全部人为生成。源数据原有的结局缺失不属于本次注入。脚本没有使用被删除的源数据基线真值来补回缺失；单位校正使用的是模拟生成的单位谱系。

本目录 CSV 被 `.gitignore` 排除。包的软件许可不能自动替代所有上游数据的再分发条款；上传前不应把患者级表直接附入仓库。原数据和模拟数据不适合用于识别个人或提供临床建议。
