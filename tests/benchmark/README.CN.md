# 基准测试：外部数据复现与嵌合体阳性对照

[English version](README.EN.md)

目的：让 MarkerFinder 有能力用**外部**数据发现自身的偏差（G6）。本目录不会自动
下载任何内容——GB 量级的下载与生信工具都由用户显式触发。

## 目录内容

| 路径 | 作用 |
|---|---|
| `expected/accessions.yaml` | GTDB ar53 目水平集合：30 个基因组（25 个目内 + 5 个外群），依 |
| `expected/chimera_positives.yaml` | 植入型嵌合体真值模板（位置已知的 HGT 阳性） |
| `expected/taxonomy_mustpass.yaml` | “must-pass” 分类学关系骨架 |
| `fetch_datasets.sh` | 从 NCBI 下载该 accession 集合到 `downloads/`（不入库） |
| `make_chimeras.py` | 构造嵌合体 FASTA：把基因组 A 的标记 X 替换为远缘基因组 B 的同源序列 |
| `run_benchmark.py` | 运行流水线、读取决策卡、计算指标 |
| `metrics.py` | precision / recall / F1 + 块自助法置信区间（块 = 基因组） |
| `provenance.py` | ：DOI / 许可证 / 获取日期 / 校验和记录；不完整时 `run_benchmark` 拒绝运行 |
| `ablation.py` | ****：反向消融——只允许两个排序在**唯一一个**因子上不同，否则拒绝 |
| `perf_probe.py` | 运行时长/规模探针，作用于随包的骨架数据 |
| `test_metrics_selfcheck.py` | pytest：指标自检，含 must-fail 对照 |
| `test_provenance_contract.py` | pytest： 门禁是双向的（随包文件必须 FAIL，填好的必须 PASS） |
| `test_ablation_contract.py` | pytest： 拒绝规则 + tau 对照 + 驱动接线 |

## 如何运行（由用户触发）

```bash
# 1. 下载基因组（GB 量级，需要网络 + 约 2 GB 磁盘）
bash tests/benchmark/fetch_datasets.sh

# 2. 构造植入型嵌合体
python tests/benchmark/make_chimeras.py --downloads tests/benchmark/downloads

# 3. 运行基准（需要 mafft/trimal/FastTree/ASTRAL 在 PATH 中）
python tests/benchmark/run_benchmark.py --downloads tests/benchmark/downloads

# 4. 反向消融——两次已完成的运行，且只声明一个改变因子。
#    该阶段不需要比对工具：它只读两个排序 TSV 并比较。
python tests/benchmark/ablation.py \
    --run-a run_with_filter/Phase5_reports/markerfinder.marker_summary.tsv \
    --run-b run_without_filter/Phase5_reports/markerfinder.marker_summary.tsv \
    --factor hgt_filter=on:off --out-dir ablation_out

# 或者通过驱动脚本，它会在探测工具之前打印同样的结论：
python tests/benchmark/run_benchmark.py --ablation-a ... --ablation-b ... \
    --ablation-factor hgt_filter=on:off
```

值得了解的 规则：省略 `--factor` 会被拒绝（没有指明原因的排序差异不构成
消融）；传两个 `--factor` 同样被拒绝——差异将无法归因。声明的设置完全相同也会
被拒绝，因为在本代码库中这种情形历史上意味着该开关是死的。分歧以“平局校正后的
Kendall tau-b”以及逐标记的秩位移报告在 `<prefix>.ablation.tsv` 中；得分为 `NA`
的标记会被剔除并在 stderr 上点名，而不会被当作零分参与排序。

如果缺少 mafft/trimal/FastTree/ASTRAL，`run_benchmark` 会逐阶段打印
`NOT EXECUTED` 并以非零码退出——一个在缺少工具时仍悄悄“通过”的基准，正是
中我们绝不可复现的失败模式。

## 来源 / 许可证 / 引用

* 基因组：`expected/accessions.yaml` 中列出的 NCBI RefSeq accession
  （选取规则见：一个 GTDB ar53 目，目内 25 个代表基因组 + 目外 5 个外群代表）。
* 分类学真值：GTDB r53（Parks 等 2021，doi:10.1038/s41587-021-00960-6）。
* 下载得到的数据**不进入版本库**；`downloads/` 已被 git 忽略，
  `fetch_datasets.sh` 会把每个文件的 SHA-256 记入 `downloads/manifest.sha256`。
