# 基准测试：外部数据复现与嵌合体阳性对照

[English version](README.EN.md)

目的：让 MarkerFinder 有能力用**外部**数据发现自身的偏差。本目录不会自动
下载任何内容——GB 量级的下载与生信工具都必须由用户显式触发。

## 目录内容

| 路径 | 作用 |
|---|---|
| `expected/accessions.yaml` | GTDB ar53 目阶元集合：30 个基因组（25 个目内 + 5 个外群），依据见[来源 / 许可证 / 引用](#来源--许可证--引用) |
| `expected/chimera_positives.yaml` | 植入型嵌合体真值模板（位置已知的 HGT 阳性） |
| `expected/taxonomy_mustpass.yaml` | “must-pass” 分类学关系骨架 |
| `fetch_datasets.sh` | 从 NCBI 下载该 accession 集合到 `downloads/`（不入库） |
| `make_chimeras.py` | 嵌合体真值门禁：校验 `expected/chimera_positives.yaml` 与 `downloads/manifest.sha256`，使替换步骤不可能静默空转；逐 marker 替换本身（把基因组 A 的 marker X 换成远缘基因组 B 的直系同源序列）**未**实现 |
| `run_benchmark.py` | 探测外部工具、下载清单与溯源记录；它所门禁的“流水线 → 决策卡 → 指标”阶段**未**实现，因此即便所有门禁都通过，最终仍以 `NOT EXECUTED` 结束 |
| `metrics.py` | precision / recall / F1 + 块自助法置信区间（块 = 基因组） |
| `provenance.py` | [`expected/accessions.yaml`](expected/accessions.yaml) 的 `provenance:` 段：DOI / 许可证 / 获取日期 / 校验和记录；不完整时 `run_benchmark` 拒绝运行 |
| `ablation.py` | **[如何运行](#如何运行由用户触发)第 4 步的反向消融**：比较两个排序，且只允许它们在**唯一一个**因子上不同，否则拒绝 |
| `perf_probe.py` | 性能预算探针：在合成数据（60 个 marker x 60 个类群，种子 `20260920`）上为新增的纯 Python 路径计时，并与 60 秒上限比较；它不读随包数据、不调用外部工具，而把“占运行总时长 5%”的分母留给用户自己的日志 |
| `test_metrics_selfcheck.py` | pytest：指标自检，含 must-fail 对照 |
| `test_provenance_contract.py` | pytest：[`provenance.py`](provenance.py) 门禁是双向的（随包文件必须 FAIL，填好的必须 PASS） |
| `test_ablation_contract.py` | pytest：[`ablation.py`](ablation.py) 拒绝规则 + tau 对照 + 驱动接线 |

## 如何运行（由用户触发）

```bash
# 1. 下载基因组（GB 量级，需要网络 + 至少 5 GB 磁盘余量；低于 5 GB 时脚本拒绝运行并以 2 退出）
bash tests/benchmark/fetch_datasets.sh

# 2. 构造植入型嵌合体
python tests/benchmark/make_chimeras.py --downloads tests/benchmark/downloads

# 3. 运行基准测试（需要 mafft/trimal/FastTree/fasttree/astral 在 PATH 中；iqtree3 不在这个必需列表里）
python tests/benchmark/run_benchmark.py --downloads tests/benchmark/downloads

# 4. 反向消融——两次已完成的运行，且只声明一个改变因子。
#    该阶段不需要比对工具：它只读两个排序 TSV 并比较。
python tests/benchmark/ablation.py \
    --run-a run_with_filter/Phase5_reports/markerfinder.marker_summary.tsv \
    --run-b run_without_filter/Phase5_reports/markerfinder.marker_summary.tsv \
    --score-column marker_quality_score \
    --factor hgt_filter=on:off --out-dir ablation_out

# 或者通过驱动脚本，它会在探测工具之前打印同样的结论：
python tests/benchmark/run_benchmark.py --ablation-a ... --ablation-b ... \
    --ablation-factor hgt_filter=on:off
```

[`ablation.py`](ablation.py) 值得了解的规则：省略 `--factor` 会被拒绝（没有指明原因的排序差异不构成
消融）；再传一个指明不同因子的 `--factor` 同样被拒绝——差异将无法归因。声明的
设置完全相同也会被拒绝，因为在本代码库中这种情形历史上意味着该开关是死的——
`ablation.py` 点出了三个曾这样死掉的开关：`--scan-stability-min`、
`--min-informative-sites` 与 `composition` 步骤。分歧以“平局校正后的 Kendall tau-b”
以及逐 marker 的秩位移报告在 `<prefix>.ablation.tsv` 中；得分为 `NA` 的
marker 会被剔除并在 stderr 上点名，而不会被当作零分参与排序。`--score-column`
默认为 `overall_score`，而 `markerfinder.marker_summary.tsv` 并没有这一列，所以
针对真实产物必须显式给出；上文的驱动脚本形式不会转发该选项，因此对真实的
`marker_summary.tsv` 会拒绝运行。

如果缺少 mafft/trimal/FastTree/fasttree/astral 中的任何一个，`run_benchmark` 会逐阶段打印
`NOT EXECUTED` 并以非零退出码退出——一个在缺少工具时仍悄悄“通过”的基准测试，正是
我们绝不可复现的失败模式。

## 来源 / 许可证 / 引用

* 基因组：`expected/accessions.yaml` 中列出的 NCBI RefSeq accession
  （选取规则见 [`expected/accessions.yaml`](expected/accessions.yaml)：一个 GTDB ar53 目，目内 25 个代表基因组 + 目外 5 个外群代表）。
* 前置条件：随包的 `expected/accessions.yaml` 目前是占位模板——`order` 仍为
  `REPLACE_WITH_ORDER_NAME`，且 `provenance:` 下五项（`source_registry`、`doi`、
  `license`、`retrieved_utc`、`checksum_sha256`）均为空。选定具体的 GTDB ar53 目、
  填入 accession 与这五项之前，`fetch_datasets.sh` 会直接拒绝下载并以 2 退出，
  `run_benchmark` 同样会因溯源记录不完整而拒绝运行。而即便所有门禁都通过，也仍
  得不到任何数值：替换与计分阶段并不存在，`run_benchmark` 最终仍会打印
  `NOT EXECUTED` 并以 2 退出（见上表对应两行）。
* 分类学真值：GTDB r53（Parks 等 2021，doi:10.1038/s41587-021-00960-6）。
* 下载得到的数据**不进入版本库**；`downloads/` 已被 git 忽略，
  `fetch_datasets.sh` 会把每个文件的 SHA-256 记入 `downloads/manifest.sha256`
  （该脚本只调用 `shasum -a 256`，没有 `sha256sum` 回退；Linux 上的 `shasum`
  是 Perl 附带的包装脚本——装了 `perl` 就有，最小化镜像里则没有——因此在缺少它的
  机器上，脚本自身的 `set -euo pipefail` 会在写清单这一步中止运行。）
