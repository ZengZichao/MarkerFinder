# MarkerFinder 测试方案与验证方法论

适用版本：MarkerFinder 0.1.0
[English version](TESTING.EN.md) · 配套：[TEST-REPORT.CN.md](TEST-REPORT.CN.md)

---

## 1. 测什么，以及为什么是这种形状

MarkerFinder 公开发布的论断有三类，每类需要不同的instrument（检验手段）：

| 论断类别 | 例子 | 检验手段 |
|---|---|---|
| 数值层 | 两个拓扑的 `rf_distance` 落在 0–1 | 单元测试，`tests/unit/` |
| 结构层 | `run_config.json` 记录了数据库版本 | 集成测试（mock 外部工具），`tests/integration/` |
| 结论层 | “HGT 筛选能去除系统发育上误导性的标记” | **针对阳性对照的真实端到端运行**，`validation/` |

结论层论断无法靠 mock 建立：一旦 mock 掉 MAFFT，比对结果就是 mock 想返回什么就
是什么。促成本项目自检层的正是撤稿事件：*Science* 2025 年的海绵/栉水母研究
（Steenwyk & King 2025, *Science* 390:751-756,
doi:10.1126/science.adw9456，2026-02-05 撤稿）正是在这个
缺口上失败的——所有内部检查都通过，而拓扑在生物学上荒谬。因此测试被拆成：

```
tests/       991 个快速用例，不需要外部工具，不下载数据
validation/  153 个在真实基因组上用真实工具跑的端到端用例
```

两层都随发布包提供；`validation/` 还额外提供其**数据**、**溯源记录**以及发布那次
运行的**归档结果**。

## 2. 能力面：什么叫“覆盖全部声称支持的功能”

“覆盖所有功能与全流程”只有在被声称的东西存在一个权威清单时才是可检验的。这个
权威清单就是代码：

* `markerfinder.cli.parser._build_parser` → 68 个命令行选项（含 4 个已文档化
  的废弃别名与 1 个已文档化的空操作开关）；
* `markerfinder.cli.parser._SUBCOMMAND_SPEC` → 4 个可续跑的分步命令；
* `markerfinder.cli.constants` → 7 个互不相同的退出码；
* `validation/capabilities.py` 中的产物清单 → 19 类输出产物
  （`Phase5_reports/*`、`Phase5_evidence/*`、`Phase4_trees/*`、
  `Phase4_alignments/*`、`Phase4_intermediate/*`、`.markerfinder/*`）；
* 9 条流程/数据性质：占居率选择、HGT 分级、状态持久化、确定性、溯源记录、
  不可测不得写成占位值、失败必须大声、数据溯源完整性、清单校验和。（分步与
  一次性运行的等价性由相关选项的用例断言，不单独列为一条性质。）

覆盖率不是一份手工维护的表。每个用例自己声明覆盖了什么：

```python
@pytest.mark.capability("hgt_scan", "product:Phase5_reports.threshold_scan.tsv")
def test_threshold_scan_writes_its_own_artifact(...):
```

`conftest.pytest_collection_modifyitems` 在用例收集阶段把这些标记导出到
`validation/.work/capability_snapshot/*.json`，随后
`test_v90_capability_coverage.py` 从两个方向把关：

1. 能力面上有条目而无用例 → **失败**；
2. 用例标记了能力面之外的条目 → **失败**（选项已被删除而它的测试标记还留着）。

矩阵会导出到 `validation/results/capability_matrix.tsv`，逐条列出覆盖每个能力的
用例 id，因此报告里的覆盖率结论与测试把关的是同一份数据。若只收集了部分用例
（例如 `pytest -k something`），该检查会显式 skip，而不是给出一个被削弱的通过。

选项的*角色*同样被记录：废弃别名的覆盖要求是“仍然可用”**且**“自我声明已废弃”；
`--no-hgt-adaptive-thresholds` 作为已文档化的空操作来覆盖，即断言它与默认行为
一致。

## 3. 测试数据

### 3.1 溯源

所有序列数据都是公开且可再生的，没有任何数值是手工录入的。

| 项目 | 来源 | 记录位置 |
|---|---|---|
| 12 个蛋白质组 | NCBI RefSeq，`datasets download genome accession … --include protein` ——由 `01_fetch_genomes.py` 下载；`genomes/` 不入库，由 `MANIFEST.sha256` 校验 | `validation/data/genomes/` + `PROVENANCE.tsv` |
| 物种名 + 7 个分类阶元 | NCBI `datasets summary genome accession` → `datasets summary taxonomy taxon <tax_id>` | `PROVENANCE.tsv` 的 `domain…species`、`tax_id` 列 |
| 组装统计（contig 数、长度、GC、蛋白数） | 同上 NCBI 报告 | `PROVENANCE.tsv` |
| 完整度 / 污染率 | NCBI 为每个组装记录的 `checkm_info` | `PROVENANCE.tsv`，并导出为 `data/checkm/*.tsv` |
| 校验和 | 每文件 SHA-256 | `PROVENANCE.tsv` + `data/MANIFEST.sha256` |
| 标记 HMM profile | `db/gtdb_markers/bac120`，由 `scripts/fetch_marker_db.py` 在本地拼装（第三方模型，不入库） | 原样复制到 `data/hmms/`（也是本地副本，不入库） |

三个构建脚本（`01_fetch_genomes.py`、`02_build_marker_sets.py`、
`03_build_fixtures.py`）可从 NCBI 重建整个数据包；它们幂等且支持 `--force`，
`validation/data/.hit_cache/` 保存其依赖的 HMM 扫描缓存。

### 3.2 数据层抓出的一个错误

本数据集的早期版本手工填写 taxonomy 表，把 `GCF_000008105.1` 标成
`o__Bacillales; f__Bacillaceae; g__Bacillus`。NCBI 对该组装的记录是
**Salmonella enterica** serovar Choleraesuis str. SC-B67——Pseudomonadota 门、
Enterobacterales 目。由于 taxonomy 表直接驱动 MAD 定根 + 单系性筛选，错误的阶元
不是笔误而是失效的实验：它决定了“哪个支系应当是单系的”。

现在有两道对策：

* lineage 由 `01_fetch_genomes.py` 从 NCBI 获取，没有 `tax_id` 就无法产表；
* `test_shipped_genomes_are_the_organisms_the_provenance_claims` 会读取每个蛋白质
  组自身的 FASTA 头（NCBI 为每个蛋白写入 `[organism]`），并要求记录的属名是其中
  的多数物种。

### 3.3 标记输入是“布局重构”，并且如实标注

`--marker-mode gtdb_tk` 消费 GTDB-Tk 的逐 marker 未比对 FASTA。GTDB-Tk 需要约
100 GB 的参考树数据包，无法随测试包发布，因此随包的 marker 文件是**对该布局的
重构**：对 `db/gtdb_markers/bac120` 中每个 profile，取每个基因组 bitscore ≥ 20
的最佳 hmmsearch 命中；一 marker 一文件，文件名主干即 marker id，tip 名为组装
accession。被保留的恰是流水线真正读取的东西（文件布局、marker id、未比对蛋白
序列、逐 marker 占居率）；**未**复现的是 GTDB-Tk 那套结合分类学的最佳命中选择，
因此这些文件不能被当作 GTDB-Tk 的发布输出。`OCCUPANCY_ALL.tsv` 记录了全部 121
个被扫描 profile 的占居率；每个 marker 目录另有本目录的 `OCCUPANCY.tsv`。

两个派生集合的存在有明确理由：

* `markers/<set>`（20 个 marker，占居率跨 0.875–1.0）与
  `markers/<set>_core`（8 个满占居 marker）——在“一切都满”的数据上无法证明阈值
  会拒绝；core 集合则保证运行够快；
* `markers/degenerate/`——从真实 marker 中截取出的 2-tip 与 1-tip 文件，并附
  README 说明其来历，用于“低于可建树最小信息量”的健壮性用例。

### 3.4 阳性对照

`markers/chimera_bench12/` 与干净集合仅有三个文件不同。每一处都把
*Bacillus anthracis* Ames（Caryophanales 目）的序列替换为 *Lactococcus lactis*
Il1403（Lactobacillales 目）的同源序列——一次位置已知的跨目转移，从被记录的候选
池中用 `random.Random(20260929)` 抽取。其余字节完全一致，因此干净运行就是这同
12 个基因组的成对阴性对照。

检测在**排除**判定处计分（`hgt_risk_level == level_3`），因为 Level 2 是流水线
仍然保留的“可疑”档：把它算作命中会虚增表观灵敏度。报告两个量：

* **差示判定**——在嵌合运行中被排除、而在干净运行中未被排除的 marker，从而背景排
  除率是被测出的，不是被假设的；
* **方向性**——每个植入 marker 在嵌合运行中的风险必须高于它在干净运行中的风险，
  且任何未植入的 marker 都不得变化。

## 4. 端到端用例的方法学

1. **只断言可观测结果。** 用例只对产物、列、退出码或实测差异下断言。“函数被调用
   过”不能作为行为论断的证据。
2. **优先差示断言。** 多数用例把流水线跑两遍、只改一个因子再比较
   （`--max-markers 2` 对 `8`，阈值 0.0 对 1.01，下限 0 对 100000）。绝对期望值会
   把数据集悄悄编码进去，差异才能隔离出旋钮本身。
3. **不改变任何东西的旋钮即失败。**
   `test_l2l3_boundary_decides_exclusion`、
   `test_occupancy_floor_rejects_markers_below_it`、
   `test_informative_site_floor_marks_markers_inconclusive_and_says_so` 等都为捕捉
   “读取了、记录了、然后丢掉”的选项而存在——这正是本代码库被系统性审计的缺陷类。
4. **失败路径是一等公民。** `test_v16_robustness_edges.py` 要求：非零退出码、不出
   现 traceback、消息点名原因、且不产出空的占位产物。
5. **受支持范围要被断言，而非假设。** 每个 marker 少于 4 个 tip 时无法推断基因树；
   用例把这一拒绝钉住（退出码 ≠ 0 且给出解释），并要求不得发布 0 字节的树文件。
6. **报告引用的数值由测试自己写入** `validation/results/metrics.json`，因此文档无
   法悄悄偏离实测。
7. **skip 必须有理由，且永不等于通过。** 被门禁挡下的判据会带着拒绝原文 skip；能
   面矩阵仍要求该选项被某个用例覆盖。

## 5. 环境

必需：Python 3.10 及以上，不设上界。声明的运行期依赖（biopython、pyyaml、ete3、
tomli），以及外部工具 `hmmsearch`、`mafft`、`trimal`、`FastTree`/`fasttree`、`astral`；
若使用 ML 基因树还需 `iqtree3`，若不跳过质量评估则需 `checkm`。

这个区间此前止于 3.12：ete3 在包导入时会 `import cgi`，而 CPython 3.13 已移除
`cgi`（PEP 594），导致 ete3 —— 以及所有 MAD / 单系性测量 —— 在更新的解释器上不可用。
`markerfinder._cgi_compat` 现在只在标准库 `cgi` 确实缺失时安装最小替身，因此 ete3 在
3.13+ 可以导入，上界也就没有继续存在的理由。

套件会自行固定环境（`validation/conftest.py::_tool_env`）：把调用它的解释器所在
`bin` 前置到 `PATH`，并把 `TMPDIR`/`MAFFT_TMPDIR` 指向 `validation/.work` 内。第二
条不是洁癖：MAFFT 用 `mktemp -dt` 建立暂存目录，失败后并不停止，而是带着空路径继续
跑，产不出比对，流水线于是报告“没有 marker 可比对”——一个不可写的 `/tmp` 就这样
变成了错误的科学结论。

发布那次运行的实测环境（解释器、依赖版本、各工具自报的版本、主机、CPU 数）归档在
`validation/results/environment.txt`。

## 6. 如何运行

```bash
# 快速层
pytest tests                      # 约 35 秒，不需要外部工具
pytest tests/benchmark            # 外部基准驱动的契约测试

# 验收层（真实流水线、真实数据）
python validation/run_validation.py -n 8            # 不含慢速用例
python validation/run_validation.py --all -n 8      # 含检测性能指标
python validation/run_validation.py --prepare       # 先从 NCBI 重建数据

# 单个功能域
pytest validation/cases/test_v06_hgt_screening.py -k composition -n 4
```

其他解释器通过在同一套套件下重复运行来检验：3.10、3.11、3.12、3.13 是 CI 矩阵，
声明区间本身也在测试之内——`tests/unit/test_supported_range_is_earned.py` 断言
manifest 不带上界，`tests/unit/test_self_test_environment.py` 断言 `--check` 与之一致。

## 7. 结果、发现与当前边界

发布那次运行的实测结果——通过/失败计数、耗时、覆盖率数字、检测指标，以及本轮重构
中套件暴露出的每一项产物/结构缺陷——见
**[docs/TEST-REPORT.CN.md](TEST-REPORT.CN.md)**（英文：
[docs/TEST-REPORT.EN.md](docs/TEST-REPORT.EN.md)）。

已声明的边界（写出来，而不是藏着）：

* marker 输入是 GTDB-Tk 的**布局**重构（§3.3），因此这里的端到端结果并不度量
  GTDB-Tk 自身的 marker 选择。
* 功能类别映射的文本取自各 profile 的 `# DESCRIPTION` 行：它检验的是
  `--cog-category-map` 的通路，不是 COG 生物学。
* 检测性能建立在 12 基因组上的 3 个植入嵌合体之上。它是带已知答案的阳性对照并有
  实测背景率，但绝不是具有流行病学代表性的 HGT 样本；更强的估计要走
  `tests/benchmark/`（30 基因组、GTDB r53 真值），因其下载体量而由用户触发。
* 20-marker 的占居率跨度下界只到 0.875：自带库是保守标记库，在这组亲缘较近的基因
  组上不存在极稀疏的 marker，因此此类输入只能靠截取模拟，并已如实标注。
