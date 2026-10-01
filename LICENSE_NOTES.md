# Licensing and attribution / 许可与署名

Original code and documentation: [MIT License](LICENSE), © 2026 Xuchi Cheng.

本仓库原创代码与文档采用 [MIT 许可证](LICENSE)，版权所有者为 Xuchi Cheng，年份为 2026。

## Data and third-party materials / 数据与第三方材料

- NHEFS teaching data are obtained through `causaldata::nhefs`. Attribution:
  Hernán MA and Robins JM, *Causal Inference: What If*;
  [book and data](https://miguelhernan.org/whatifbook),
  [causaldata](https://github.com/NickCH-K/causaldata).
- The installed `causaldata` version declares `MIT + file LICENSE` for the
  package. This does not establish redistribution terms for all upstream
  data. This repository's MIT license does not relicense upstream data,
  dependencies or embedded third-party report resources.
- Person-level source and simulated CSVs, and local RDS objects, are excluded
  from Git. Scripts regenerate the teaching inputs locally.
- Artificial source labels, extraction metadata, missingness and data damage
  are simulation features, not newly collected EHR or registry data.

NHEFS 教学数据通过 `causaldata::nhefs` 读取，保留上述教材与软件包来源署名。
软件包许可与上游数据再分发条件分别处理；本仓库 MIT 许可证不改变上游数据、
依赖程序包或报告内嵌第三方资源的许可。逐人 CSV 和本地 RDS 不随仓库上传，
通过脚本在本地再生成。人工来源标签及质量降级机制不代表新采集的临床数据。
