# MarkerFinder 验证套件

**基于真实数据的全功能、全流程验证。**
[English document](README.EN.md)

本目录是随 MarkerFinder 一起发布的验收证据。它用真实基因组和真实 HMM profile
运行真实流水线——这里没有任何 mock，也不是单元测试。快速的、基于 mock 的测试层
位于 [`../tests/`](../tests)。

> 本仓库的文档以英文版本为主，中文版本为对照译本；若两者出现不一致，以英文版本
> 为准。

---

## 本套件确立的结论

| | 被确立的论断 | 位置 |
|---|---|---|
| 1 | 命令行公开的每一项能力都有对应用例 | `cases/test_v90_capability_coverage.py`、`results/capability_matrix.tsv` |
| 2 | 全流程既可分步执行也可一次性执行，且结果一致 | `cases/test_v12_stepwise_workflow.py` |
| 3 | 每一个已文档化的输出文件都存在，且列与信息一致 | `cases/test_v14_product_audit.py` |
| 4 | 检测性能是针对"已知答案的阳性对照"实测得到的 | `cases/test_v17_benchmark_detection.py` |
| 5 | 相同输入两次运行结果一致 | `cases/test_v15_reproducibility.py` |
| 6 | 失败路径以正确的退出码大声失败，且不出现 traceback | `cases/test_v16_robustness_edges.py` |
| 7 | 测试数据本身与其声明一致 | `cases/test_v01_environment_and_selfcheck.py` |

---

## 目录结构

```
validation/
├── README.EN.md / README.CN.md   本文件（以英文版为主）
├── capabilities.py               能力清单，直接从代码读取
├── conftest.py                   fixture：运行器、输入、每用例独立输出
├── run_validation.py             一键执行入口
├── cases/                        test_vNN_*.py，按功能域编号
├── data/                         随包夹具 + 需下载的输入数据（见下文）
│   ├── PROVENANCE.tsv            每个组装的物种、tax_id、 lineage、校验和
│   ├── MANIFEST.sha256           本目录所有文件的 SHA-256 清单
│   ├── genomes/                  12 个 NCBI RefSeq 蛋白质组（.faa）——需下载，不入库
│   ├── sets/                     各命名输入集合包含哪些 accession
│   ├── markers/                  逐 marker FASTA、占居率表、嵌合体对照
│   ├── hmms/                     hmm 模式使用的 TIGRFAM/Pfam profile
│   ├── taxonomy/                 Format A / Format B 表与畸形行对照
│   ├── trees/                    参照物种树与非法参照树
│   ├── checkm/                   取自 NCBI 自身记录的质量评估表
│   ├── cog/                      marker → 功能类别映射
│   ├── mustpass/                 分类学 must-pass 基线（通过与失败两向）
│   ├── configs/                  YAML / TOML / JSON 配置文件
│   ├── sequences/                蛋白与核酸附属文件
│   └── variants/mag_named/       用于 MAG 路径的标识符命名变体
├── scripts/                      数据如何从公开来源构建
│   ├── 01_fetch_genomes.py       NCBI 下载 + 权威 lineage 获取
│   ├── 02_build_marker_sets.py   逐 marker FASTA、占居率、植入型嵌合体
│   └── 03_build_fixtures.py      表、树、质量文件、配置、清单
└── results/                      发布那次运行的归档产物
    ├── metrics.json              测试报告中引用的全部数值
    ├── capability_matrix.tsv     能力 → 覆盖用例
    ├── chimera_metrics.json      检测的 precision / recall / F1
    ├── junit.xml                 逐用例通过失败记录
    └── environment.txt           实测的解释器、依赖包与工具版本
```

`.work/` 是临时区：materialize 出来的输入目录、各用例输出与临时目录。它可随时
删除，且不进入版本库。

---

## 如何运行

前置条件：MarkerFinder 环境（Python 3.10–3.12，外加 `hmmsearch`、`mafft`、
`trimal`、`FastTree`/`fasttree`、`astral`、`iqtree3`），并执行
`pip install -e ".[dev]"`。

```bash
# 除慢速 12 基因组检测用例外的全部用例，8 个并行 worker
python validation/run_validation.py -n 8

# 包含慢速用例（在 12 个基因组上实测检测指标）
python validation/run_validation.py --all -n 8

# 先从 NCBI 重建测试数据，再运行
python validation/run_validation.py --prepare
```

等价的直接 pytest 调用：

```bash
pytest validation/cases -m "not slow" -n 8
pytest validation/cases -m slow            # 检测指标，约 10 分钟
pytest validation/cases/test_v06_hgt_screening.py -k composition
```

套件会自行固定运行环境：把解释器所在的 `bin` 前置到 `PATH`，并把
`TMPDIR`/`MAFFT_TMPDIR` 指向 `validation/.work` 内。否则 `/tmp` 只读或写满时，
MAFFT 会静默失败并表现为"没有任何 marker 可比对"，把一个环境问题伪装成科学结论。

### 运行时长

在 96 核机器上以 `-n 14` 实测：非慢速套件约 15 分钟完成。单次流水线运行在默认
"4 基因组 / 4 marker"配置下约 40 秒；同一数据集换成 30 个 marker 则需要
7 分 19 秒。因此功能级用例只使用小 marker 预算，仅占居率与检测性能用例使用较大
集合。

---

## 关于测试数据的如实说明

**基因组。** 12 个完整的 NCBI RefSeq 组装，蛋白质 FASTA（`--include protein`），
覆盖 Bacillota 门的两个目，外加一个 Pseudomonadota 外群。
`data/PROVENANCE.tsv` 逐条记录：物种名、tax_id、NCBI 的七个分类阶元、组装级别、
contig 数、序列长度、GC%、蛋白数、NCBI 自身的 CheckM 完整度/污染率、下载地址、
获取日期以及随包文件的 SHA-256。

分类 lineage 是**从 NCBI 读取的，不是手打的**。本数据集的早期版本手工填写
taxonomy 表，把 *Salmonella enterica*（`GCF_000008105.1`）标成
`o__Bacillales;f__Bacillaceae;g__Bacillus`。由于 taxonomy 表直接驱动单系性筛选，
那是一个失效的实验，而不是一个笔误。现在
`test_shipped_genomes_are_the_organisms_the_provenance_claims` 会用每个蛋白质组
自身的 header 校验其记录的属名，而 `01_fetch_genomes.py` 不产生任何没有
`tax_id` 的表格。

**逐 marker FASTA（`markers/`）。** GTDB-Tk 会在 `align/marker_genes/` 下为每个
marker 写一个未比对的 FASTA，`--marker-mode gtdb_tk` 消费的正是这种布局。
GTDB-Tk 本身需要约 100 GB 的参考树数据包，无法随测试包发布，因此这些文件是
**对该布局的重构**：用 `db/gtdb_markers/bac120` 中的 TIGRFAM/Pfam（由 `scripts/fetch_marker_db.py` 拼装）
profile 对每个基因组做 hmmsearch，取 bitscore ≥ 20 的最佳命中；一 marker 一文件，
文件名主干即 marker id，tip 名为组装 accession。被保留的是流水线真正读取的一切；
未被复现的是 GTDB-Tk 自身那套结合分类学的最佳命中筛选。
`markers/OCCUPANCY_ALL.tsv` 记录了所有被扫描 profile 的占居率，每个目录另有
自己的 `OCCUPANCY.tsv`。

**植入型嵌合体（`markers/chimera_bench12`）。** 在 3 个 marker 中，把
*Bacillus anthracis* Ames 的序列替换为 *Lactococcus lactis* Il1403 的同源序列，
即一次位置已知的跨目转移；候选池与抽取规则（`random.Random(20260929)`）都记录在
`markers/chimera_ground_truth_bench12.json` 中。目录内其余字节与干净集合完全一致，
因此干净运行同时充当这 12 个基因组的阴性对照。

**质量数值（`checkm/`）。** 完整度与污染率取自 NCBI 每个组装自己的 `checkm_info`
数值，不是编造的。NCBI 未记录数值的组装会被省略，且省略原因写入文件本身。

**功能类别（`cog/`）。** 类别文本取自该 marker 的 HMM profile 的 `# DESCRIPTION`
行。它验证的是 `--cog-category-map` 的通路是否真的接通，并不对 COG 生物学做任何
断言，文件里也这样写明。

**MAG 命名变体（`variants/mag_named/`）。** 把 quad4 的 4 个基因组改名为
`mag_01..mag_04`，并同步改写 marker 头部与 taxonomy 表。序列完全未变——分类规则
依据文件名生效，所以这是一个标签实验，并且被如实记录为标签实验。

---

## 用例遵循的约定

* 每个用例都有独立输出目录，因此套件可安全地以 `-n` 并行执行。
* 用例只断言**可观测**的结果：某个产物、某列、某个退出码，或"只改一个旋钮"的
  两次运行之间的差异。若某个旋钮不改变任何东西，就被记为失效对照，而不是通过。
* 优先使用差异断言（A/B 两次运行，只差一个因子），因为绝对阈值会把数据集本身
  悄悄编码进期望值。
* skip 必须写明理由；skip 永远不会被当作 pass。
* 报告引用的数值由各用例自己写入 `results/metrics.json`，因此文档里出现陈旧数字
  会表现为机读记录中的差异，而不是无人察觉的漂移。
