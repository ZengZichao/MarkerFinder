# MarkerFinder

**自适应 HGT 感知的原核生物系统发育基因组学流水线**

[![License: MIT](https://img.shields.io/badge/License-MIT-blue.svg)](LICENSE)
[![Version](https://img.shields.io/badge/version-0.1.0-green.svg)](https://github.com/ZengZichao/MarkerFinder/releases)
[![Python 3.10～3.12](https://img.shields.io/badge/python-3.10--3.12-blue.svg)](https://www.python.org/downloads/)
[![Tests](https://img.shields.io/badge/tests-988%20collected-brightgreen.svg)](#测试)

[English](README.EN.md) | [详细手册（EN）](MANUAL.EN.md) | [详细手册（CN）](MANUAL.CN.md)

---

## 概述

MarkerFinder 是一条面向原核生物的系统发育基因组学流水线，用于标记基因选择与物种树推断。它针对现有方法的三类常见问题：

1. **标记基因集僵化** — 根据输入数据集动态选择最优标记组合
2. **HGT 污染未处理** — 用基于系统发育的水平基因转移筛查（Phylogenetic），剔除会误导系统发育分析的标记
3. **MAG 支持不足** — 用质量感知的自适应参数，稳定处理不完整的宏基因组组装基因组

**生物学领域：** 原核生物系统学、系统发育基因组学、宏基因组组装基因组（MAG）分类、微生物分类学。

---

## 流水线架构

```
蛋白序列输入 (.faa)
        │
        ▼
┌─────────────────────────────────────────────────────┐
│  Phase 0: 质量感知预处理                              │
│  CheckM 质量评估 → 质量分层                           │
│  → 自适应参数计算                                      │
└───────────────────────┬─────────────────────────────┘
                        ▼
┌─────────────────────────────────────────────────────┐
│  Phase 1: GTDB-TK 标记基因加载与筛选                    │
│  消费 GTDB-TK ar53/bac120 每标记序列 → 占有率矩阵        │
│  → 基于占有率选择标记子集                               │
└───────────────────────┬─────────────────────────────┘
                        ▼
┌─────────────────────────────────────────────────────┐
│  Phase 2: HGT 筛查（系统发育）                              │
│  系统发育步骤: MAD 定根 + 单系比例                              │
│  或 RF距离 / Quartet一致性与物种树比较                          │
│  → 等级 1（清洁）/ 等级 2（可疑）/ 等级 3（排除）       │
└───────────────────────┬─────────────────────────────┘
                        ▼
┌─────────────────────────────────────────────────────┐
│  Phase 3: 双策略系统发育推断                            │
│  模式A: MAFFT → trimAl → IQ-TREE3 (串联法)             │
│  模式B: MAFFT → trimAl → 基因树 → ASTRAL-III (合并法)   │
│  → 冲突检测 (RF距离、Quartet分析)                       │
└───────────────────────┬─────────────────────────────┘
                        ▼
┌─────────────────────────────────────────────────────┐
│  Phase 4: 报告生成                                    │
│  HTML 静态报告 + 纯文本输出                           │
│  (TSV、Newick、Nexus) + Phase5_metadata/run_config.json       │
└─────────────────────────────────────────────────────┘
```

> 注：合并法（模式 B）的基因树默认由 FastTree + WAG 构建；可通过 `--gene-tree-builder iqtree` 切换为 IQ-TREE3。

---

## 依赖要求

### Python 依赖

| 包名 | 最低版本 | 用途 |
|-----------|---------|--------------------------------------------------------|
| Biopython | >= 1.81 | 序列解析 |
| PyYAML | >= 6.0 | 配置文件解析 |
| ete3 | >= 3.1 | 系统发育树操作（Python 3.10～3.12；3.13 及以上需自行回退） |
| tomli | >= 2.0 | TOML 配置解析（Python < 3.11） |

> 注：本项目此前声明了 `pandas`、`numpy`、`scipy` 和 `rich` 四个依赖，包内却从未导入它们，
> 现已从运行时依赖中移除。`--check` 只把它们作为提示信息报告，不参与判定。

### 环境与验证（受支持解释器）

`pyproject.toml` 声明 `requires-python = ">=3.10,<3.13"`。**ete3 在 Python 3.13 及以上无法导入**：
它的 `ete3/webplugin/webapp.py` 仍会 `import cgi`，而标准库的 `cgi` 已在 3.13 中移除。
因此 RF/quartet 的 ete3 路径只在 3.10～3.12 上可用。在 3.13 及以上，`--check` 会明确报告
“版本越界”并以非零码退出，同时启用纯标准库的 split-set 回退（`method="splits-python"`）。

推荐用仓库根目录的 `environment.yml` 创建 conda 环境，它会一并安装 mafft、trimal、FastTree、
IQ-TREE 和 ASTRAL：

```bash
conda env create -f environment.yml
conda run -n markerfinder markerfinder --check
conda run -n markerfinder pytest -q
```

**使用 pip 安装时需注意**：PyPI 上的 `ete3` 未声明依赖 `six`。只执行 `pip install ete3` 时，
`import ete3` 会抛出 `ModuleNotFoundError: No module named 'six'`。请改用
`pip install ete3 six`。conda 渠道的 `ete3` 没有此问题。

**受支持范围内的实测结果**（同一份代码，在每个受支持解释器上配合 ete3 完整运行一次）：

| 解释器 | 结果 |
|---|---|
| CPython 3.10.20 / 3.11.15 / 3.12.13 + ete3 | 全量 **0 skipped / 0 failed**（用例总数以顶部徽章为准，随提交实跑刷新，此处不重复）；`--check` 49 项检查 0 FAIL、退出码 0（未获取参考数据库的源码检出为 48 项检查，那一项以 INFO 报出） |
| CPython 3.12.13 | ete3 在该版本会打印大量 `SyntaxWarning: invalid escape sequence`。这些警告与被测代码无关，但会干扰用户阅读 |
| CPython 3.14.6（**超出声明范围**） | ete3 不可导入：依赖 ete3 的模块以含 `NOT EXECUTED` 的理由显式 skip，其余用例全部通过、退出码 0；`--check` 明确报告“版本越界”并以非零码退出 |

ete3 差分测试（把纯 Python split-set 的 RF/quartet 结果与 ete3 逐树对比）**只有在
3.10～3.12 + ete3 的环境下才会真正执行**。在无 ete3 的环境里，它们以含 `NOT EXECUTED` 的理由
显式 skip，不计为通过。

声明范围的上界本身也在测试之内（`test_supported_range_is_earned.py`）。如果在 3.13 及以上，
ete3 反而变得可以导入，该测试就会失败，并要求**重新决定** `<3.13` 这个上界，
而不是让一条边界条件依赖记忆维持。

在 Python 3.14（超出声明范围）上直接运行 `pytest tests` 时，两个依赖 ete3 的模块
会在收集阶段被显式跳过并标注 `NOT EXECUTED`（而不是崩溃后一条结果都不报，退出码 2）；
3.11 与 3.12 已一并纳入实测范围。

### 外部工具

| 工具 | 最低版本 | 用途 |
|-------------|----------|----------------------------------------|
| MAFFT | >= 7.520 | 多序列比对（标记内的物种间比对） |
| trimAl | >= 1.4 | 比对修剪 |
| IQ-TREE3 | >= 3.0.0 | 最大似然建树 |
| ASTRAL-III | >= 5.7.0 | 合并法物种树 |
| FastTree2 | >= 2.1.11| 快速基因树构建（默认 FastTree + WAG） |
| DIAMOND | >= 2.1.0 | BBH 直系同源物解析（当前主流程跳过，接口预留）|
| CheckM | >= 1.2 | 基因组质量评估（可选） |

### 数据库

| 数据库 | 说明 | 路径 | 是否必需 |
|-------------|--------------------------------|-------------------|---------|
| GTDB-TK 提取的每标记 FASTA | bac120（细菌）与 ar53（古菌）保守蛋白，每个标记一个文件（文件名主干即标记 id） | `--gtdb-markers-dir` | `gtdb_tk` 模式必需 — MarkerFinder 直接消费 GTDB-TK 输出 |
| TIGRFAM/Pfam 单标记 HMM | `hmm` 模式下 MarkerFinder 自动从 `--marker-hmm-dir` 或 `--db-dir` 下的 `gtdb_markers/{ar53,bac120}` 发现。**本仓库不携带这些文件**——它们是第三方模型，请用 `python scripts/fetch_marker_db.py` 在本地拼装（见 [`db/README.md`](db/README.md)） | `--marker-hmm-dir` | `hmm` 模式必需（未指定时尝试自动发现） |

> 默认标记发现模式为 **`gtdb_tk`**。`hmm` 模式用于不依赖 GTDB-TK 预运行的场景。
> `--marker-mode` 只接受 `gtdb_tk` 与 `hmm` 两个值，不接受 `denovo`。
> HGT 筛查只基于**系统发育**（Phylogenetic）步骤，**无需任何外部大型蛋白数据库**（例如 40 GB～80 GB 的 NCBI nr），整套流程自包含运行。
> 数据库版本由 SHA256 哈希校验锁定。运行时主要校验 `gtdb_markers` 目录（`gtdb_tk` 模式）或 `marker_hmm_dir`（`hmm` 模式），版本信息记录在 `Phase5_metadata/run_config.json` 中，用于保证可复现性。

---

## 安装

### 方式一：Conda（推荐）

```bash
# 创建环境
conda create -n markerfinder python=3.10 -y
conda activate markerfinder
# 或直接基于 environment.yml 创建：conda env create -f environment.yml

# 安装外部工具
conda install -c bioconda -c conda-forge \
    mafft trimal iqtree astral fasttree \
    diamond checkm-genome ete3

# 1. 先用 GTDB-TK 对基因组集提取保守蛋白序列与物种树 (一次运行, MarkerFinder 不调用 GTDB-TK)
#    https://ecogenomics.github.io/GTDBTk/

# 安装 MarkerFinder
pip install markerfinder

# 或从源码安装
git clone https://github.com/ZengZichao/MarkerFinder.git
cd markerfinder
pip install -e ".[dev]"
```

### 方式二：从源码安装

```bash
git clone https://github.com/ZengZichao/MarkerFinder.git
cd markerfinder
pip install -e ".[dev]"
```

### 验证安装

```bash
# 检查 CLI
markerfinder --version
# 预期输出: markerfinder 0.1.0
# 版本号是包内的常量，因此解压的源码树、可编辑安装和构建出的 wheel 报的是同一个串。

# 检查外部工具
mafft --version                # MAFFT
iqtree3 --version              # IQ-TREE3
diamond version                # DIAMOND

# 运行单元测试
pytest tests/unit -q
```

---

## 快速开始

### 最小可运行示例

```bash
# 1. 准备输入目录，包含蛋白 FASTA 文件 (.faa/.fasta/.fa)
mkdir -p example/genomes
# 每个 .faa 文件应包含一个基因组的蛋白序列
# 以下命令以 gtdb_tk 模式为例，需用户先通过 GTDB-Tk 准备 --gtdb-markers-dir

# 2. 使用标准模式运行 MarkerFinder
markerfinder \
    -i example/genomes \
    -o example/output \
    -t 8 \
    --mode standard \
    --gtdb-markers-dir gtdbtk_out/align/marker_genes

# 3. 查看结果
ls example/output/
# Phase5_reports/markerfinder.report.html         — 静态 HTML 报告
# Phase5_reports/markerfinder.marker_summary.tsv  — 标记基因汇总表
# Phase4_trees/markerfinder.species_tree_concat.newick — 物种树 (Newick)
```

### MAG 分析示例

```bash
markerfinder \
    -i mags/ \
    -o output/ \
    -t 8 \
    --mode mag_adaptive \
    --gtdb-markers-dir gtdbtk_out/align/marker_genes \
    --checkm-results checkm_quality_report.tsv
```

### 跳过 CheckM（快速测试）

```bash
markerfinder \
    -i genomes/ \
    -o output/ \
    -t 8 \
    --mode standard \
    --gtdb-markers-dir gtdbtk_out/align/marker_genes \
    --skip-checkm
```

> `--skip-checkm` 会跳过 CheckM 质量评估，直接使用默认质量估计继续流程。`pipeline_summary.txt` 中 `Quality source` 会显示为 `skipped`。

### 使用 GTDB-TK 提供的保守蛋白序列

MarkerFinder 是 GTDB-TK 的下游工具，直接消费 GTDB-TK 为输入基因组集提取的
ar53（古菌）与 bac120（细菌）每标记原始序列，并把它们作为候选标记集。
**MarkerFinder 不需要再做任何 HMM 扫描或 DIAMOND 自比对**。

GTDB-TK 的 ar53/bac120 输出是**未对齐、未修剪**的原始序列。软件先用 MAFFT 逐标记比对，
再用 trimal 修剪，然后逐标记建树。HGT 系统发育步骤的判定有两种方式：

- **推荐**：提供 `--taxonomy-table`。每个标记的基因树先用 MAD（Minimal Ancestor
  Deviation）定根，再按基因树涵盖的分类学范围（scope）**自动选取下一阶元**统计单系
  比例，风险 = `1 − 单系比例`。这种方式无需物种树：域树查门、门树查纲、纲树查目、目树查科、
  科树查属、属树查种。用 `--monophyly-rank` 显式指定阶元时（默认为 `auto`），就按
  指定的那个阶元统计。只有当分类单元在 split 两侧各拥有 ≥ 2 个代表时才会计入分母；
  单系判定看的是基因树是否携带该 split（两侧任一为 clade 即算成立），因此结论不依赖定根位置。
- 或提供 `--species-tree`，把基因树与 GTDB-TK 的串联物种树比较，
  以 RF/quartet 不一致判定 HGT。

```bash
# 1) 先用 GTDB-TK 对基因组集提取保守蛋白序列与串联物种树
gtdb-tk align    --genome_dir genomes/ --out_dir gtdbtk_out/ \
                 --cpus 8 --extension gz
gtdb-tk classify --genome_dir genomes/ --out_dir gtdbtk_out/ \
                 --cpus 8

# 2a) 推荐:提供分类学表，系统发育步骤走 MAD 定根 + 单系比例（无需物种树）
markerfinder \
    -i genomes/ \
    -o output/ \
    -t 8 \
    --gtdb-markers-dir gtdbtk_out/align/marker_genes \
    --taxonomy-table gtdbtk_out/classify/gtdbtk.bac120.summary.tsv \
    --monophyly-rank auto \
    --monophyly-threshold 0.5 \
    --hgt-steps phylogenetic

# 2b) 或者：提供 GTDB-TK 串联物种树作系统发育步骤参照
markerfinder \
    -i genomes/ \
    -o output/ \
    -t 8 \
    --gtdb-markers-dir gtdbtk_out/align/marker_genes \
    --species-tree gtdbtk_out/classify/gtdbtk.bac120.classify.tree \
    --hgt-steps phylogenetic
```

> MarkerFinder 的标记集严格等于 GTDB-TK 的保守蛋白集合。GTDB-TK 的输出未对齐，
> 因此建树前先用 MAFFT 比对、再用 trimal 修剪。流程只剔除高 HGT 风险的标记，其余全部保留。
> 分析单一域时，只需提供对应的 bac120（细菌）或 ar53（古菌）目录，以及分类学表或物种树。

### 分步运行

MarkerFinder 支持把流程拆分为 4 个子命令（`scan`、`filter`、`infer`、`report`）逐步执行。
运行状态在输出目录之间自动传递：

```bash
markerfinder scan     -i genomes/ -o output/ -t 8 --gtdb-markers-dir gtdbtk_out/align/marker_genes
markerfinder filter             -o output/ -t 8
markerfinder infer              -o output/ -t 8 --gene-tree-builder fasttree
markerfinder report             -o output/
```

详见 [MANUAL.CN.md](MANUAL.CN.md)。

---

## 输入格式

输入目录必须至少包含一个**蛋白 FASTA 文件**，扩展名为 `.faa`、`.fasta` 或 `.fa`。
仅仅提供 `.fna` 核苷酸输入已不受支持。每个 `.faa` 文件对应一个基因组的预测蛋白序列。
CheckM 质量评估为可选项，可用 `--skip-checkm` 显式跳过。

```
input_dir/
├── genome_1.faa     # 基因组1的预测蛋白序列
├── genome_2.faa     # 基因组2的预测蛋白序列
├── mag_bin1.faa     # MAG蛋白序列（自动识别为MAG）
└── ...
```

**文件命名约定：**

- 文件名包含 `mag` 或 `bin` 的，自动归类为 MAG
- 文件名包含 `sag` 的，归类为单细胞扩增基因组（SAG）
- 其余文件归类为分离株基因组

**序列格式：** 标准 FASTA 格式。每个文件内，每个基因只保留一条蛋白序列，
标题行格式为 `>序列ID [可选描述]`。

---

## 核心参数

### 基本参数

| 参数 | 默认值 | 说明 |
|------|--------|------|
| `-i, --input` | *必填* | 包含 `.faa` 文件的输入目录 |
| `-o, --output` | *必填* | 输出目录 |
| `-t, --threads` | `min(8, CPU 数)` | CPU 线程数 |
| `--mode` | `standard` | 分析模式：`conservative`、`standard`、`expanded`、`mag_adaptive` |
| `--config` | — | YAML/TOML/JSON 配置文件路径（CLI 参数优先于配置文件） |

### 树与分类学参数

| 参数 | 默认值 | 说明 |
|------|--------|------|
| `--tree` | — | `--species-tree` 的别名，用于 HGT 系统发育步骤的参考物种树（RF/quartet） |
| `--sequences` | — | 参考序列文件（FASTA） |
| `--taxonomy-table` | — | 外部分类学表格（TSV/CSV） |
| `--taxonomy-format` | `table` | 分类学格式：`table`、`embedded` |
| `--taxonomy-source-priority` | `table` | 分类学来源优先级：`table`、`embedded` |
| `--taxonomy-delimiter-mode` | `reverse` | Format A 解析策略：`reverse`、`greedy`、`segment` |
| `--table-sep` | — | 强制指定表格分隔符（默认自动检测） |
| `--multi-tree-mode` | `ask` | 多棵树处理：`ask`、`split`、`first`、`last`、`random` |
| `--strip-annotations` | `false` | 去除树中的 NHX 注释 |
| `--mol-type` | `auto` | 分子类型：自动检测或 `DNA`、`RNA`、`protein` |
| `--skip-length-check` | `false` | 跳过序列长度一致性检查 |
| `--no-cross-check` | `false` | 跳过树-序列交叉验证 |
| `--ignore-malformed` | `false` | 跳过畸形行/文件而非终止（默认终止） |
| `--taxonomy-levels` | — | 以 `level:prefix` 形式（Format B 如 `kingdom:k__`，Format A 如 `kingdom:_k_`）在内置阶元之外扩展自定义阶元。该阶元会被从 `--taxonomy-table` 中解析出来并随分类学映射传递到下游，而不再以 "Unknown level prefix" 警告丢弃。注意单系筛查的阶元阶梯仍为 `domain..species`，自定义阶元只被携带、不被筛查 |

### HGT 筛查参数

| 参数 | 默认值 | 说明 |
|------|--------|------|
| `--hgt-steps` | `phylogenetic` | 逗号分隔的启用筛查步骤列表。`phylogenetic` = 基因树与参照物种树的 RF/Quartet 不一致，或在有 taxonomy 表时用 MAD 定根 + 单系比例。`composition` = 的 RCV/GC 组成偏差诊断，以**平行证据列**写入 `Phase5_reports/{prefix}.composition.tsv`，**永不**并入风险分。默认仅 `phylogenetic` |
| `--marker-mode` | `gtdb_tk` | 标记发现模式：`gtdb_tk`（默认，直接消费 GTDB-TK ar53/bac120 每标记 FASTA）或 `hmm`（基于 TIGRFAM/Pfam HMM 扫描）。只接受这两个值，无 `denovo` |
| `--gtdb-markers-dir` | — | `gtdb_tk` 模式**必填**：GTDB-TK 提取的每标记 FASTA 目录（每保守蛋白一个文件，文件名主干即标记 id）；GTDB-TK 输出为未对齐/未修剪原始序列，软件会自动 MAFFT 比对 + trimal 修剪后再建树 |
| `--marker-hmm-dir` | — | `hmm` 模式**必填**（未指定时尝试按 `--marker-db-source` 自动发现）：TIGRFAM/Pfam 单标记 HMM 目录，每个标记一个 `.HMM`/`.hmm` 文件，文件名主干即标记 id |
| `--species-tree` | （空） | `gtdb_tk` 模式可选项：GTDB-TK 串联物种树（Newick），作为 HGT 系统发育步骤的 RF/quartet 参照。**推荐**改用 `--taxonomy-table` 让系统发育步骤走 MAD 定根 + 单系比例 |
| `--taxonomy-table` | （空） | `gtdb_tk` 模式推荐：外部分类学表（TSV/CSV，首列=基因组 id），配合 MAD 定根 + 单系比例做系统发育步骤筛查 |
| `--monophyly-rank` | `auto` | `gtdb_tk` + `--taxonomy-table` 时单系比例的**测量阶元**。`auto`（默认）由基因树的分类学范围推导（取 scope 的下一阶元：域→门、门→纲、纲→目、目→科、科→属、属→种；tip 在任何阶元都不一致时取 `genus`）。显式命名阶元（`domain`/`phylum`/`class`/`order`/`family`/`genus`/`species`）则就在该阶元测量；若该阶元上没有分类单元具备信息量充分的 split（两侧各 ≥ 2 代表），回退向更高阶元，并将结果作为跨阶元比较报为 UNKNOWN，而不用作可与 `--monophyly-threshold` 比较的比例 |
| `--monophyly-threshold` | `0.5` | 单系比例下限，低于该值视为系统发育不一致（HGT 倾向），系统发育步骤风险 = `1 − 单系比例` |
| `--hgt-threshold` | `0.25` | HGT 风险阈值（等级分级） |
| `--hgt-adaptive-thresholds` | `false` | 启用自适应远缘放宽（**默认关闭**；跨科/跨目等远缘数据集因系统发育 HGT 风险趋于饱和，开启后把 `level2_max` 放宽至 0.95，以免流水线把所有标记都误判为 Level 3） |
| `--no-hgt-adaptive-thresholds` | `false` | 显式禁用自适应远缘放宽（与默认行为一致，仅供显式表达） |

### 系统发育推断参数

| 参数 | 默认值 | 说明 |
|------|--------|------|
| `--ufboot` | `1000` | IQ-TREE3 UFBOOT 重复次数 |
| `--fast-tree` | `false` | 使用 FastTree2 构建合并法基因树（遗留别名，等价于 `--gene-tree-builder fasttree`；当前默认构建器已为 fasttree） |
| `--gene-tree-builder` | `fasttree` | 合并法基因树构建器：可选 `fasttree`（FastTree2 + WAG，快，默认）或 `iqtree`（IQ-TREE3，更充分） |
| `--coalescent-mode` | `post-filter` | 合并法推断：`off`、`post-filter`、`always` |

### MAG 与输出参数

| 参数 | 默认值 | 说明 |
|------|--------|------|
| `--checkm-results` | — | 预计算的 CheckM 质量报告 |
| `--skip-checkm` | `false` | 强制跳过 CheckM，使用默认质量估计 |
| `--min-occupancy` | （空） | 显式指定标记选择的最低占用率下限（0～1），设置后覆盖 Phase 0 自适应阈值；不指定时沿用自适应分层结果 |
| `--min-marker-coverage` | （空） | **已废弃**，等价于 `--min-occupancy` 的兼容别名 |
| `--db-dir` | `./db` | 数据库目录 |
| `--report-format` | `html` | 报告格式（当前仅支持 `html`） |
| `--force` | `false` | 强制覆盖已有输出文件 |
| `--no-clobber` | `false` | 跳过已有文件，不覆盖 |
| `--save-intermediates` | `false` | 将有用中间产物复制到 `<output>/Phase4_intermediate/` |
| `--redo` | `false` | 分步模式下强制重新运行已完成的步骤 |
| `--resume` | `false` | 分步模式下自动补齐缺失的前置步骤并继续运行 |
| `--tmp-dir` | 系统临时目录下 `markerfinder-{uuid8}` | 临时文件目录（默认每趟 run 独立并自动清理；显式指定路径时不自动清理） |
| `--keep-tmp` | `false` | 保留临时目录用于调试 |
| `--log-file` | — | 日志输出文件路径 |
| `-v, --verbose` | `0` | 详细输出（可叠加：`-v` 为 INFO，`-vv` 为 DEBUG） |

### 分析模式

| 模式 | 选择策略 | 最低占有率 | 最大标记数 | 适用场景 |
|------|---------|-----------|-----------|---------|
| `conservative` | InfoMax | 0.90 | 30 | 高质量分离株，严格过滤 |
| `standard` | Greedy | 0.75 | 60 | 混合分离株数据集 |
| `expanded` | RateBalanced | 0.50 | 120 | 多样性数据集，容忍缺失 |
| `mag_adaptive` | SparseOptimized | 0.40 | 150 | MAG 为主的数据集 |

---

## 输出文件

### 静态 HTML 报告

- `Phase5_reports/markerfinder.report.html` — 自包含静态 HTML 报告，包含：
  - 执行摘要仪表盘
  - 标记基因汇总表与 HGT 风险表
  - 输出文件路径说明
  - 物种树来源说明

> 当前仅支持 HTML 静态报告，不支持 Plotly 交互式可视化或 PDF。

### 纯文本输出

输出目录按运行阶段分为子目录，子目录名称带有 `PhaseN_` 前缀以体现管线步骤：

| 文件 | 格式 | 说明 |
|------|------|------|
| `Phase5_reports/markerfinder.marker_summary.tsv` | TSV | 选定标记的占有率、质量、等级 |
| `Phase5_reports/markerfinder.hgt_evaluation.tsv` | TSV | HGT 风险评分和等级 |
| `Phase5_reports/markerfinder.pipeline_summary.txt` | TXT | 人类可读的运行摘要 |
| `Phase5_reports/markerfinder.report.html` | HTML | 静态 HTML 报告 |
| `Phase4_trees/markerfinder.species_tree_concat.newick` | Newick | 串联法物种树（主要使用 IQ-TREE3；IQ-TREE3 失败时自动回退至 FastTree2） |
| `Phase4_trees/markerfinder.species_tree_astral.newick` | Newick | 合并法物种树（仅在 ASTRAL-III 成功时输出；失败时不输出空树文件，仅记录 ERROR 日志） |
| `Phase4_trees/markerfinder.gene_trees.newick` | Newick | 各标记基因树集合（标准 multi-newick，每行一棵树），启用合并法时生成 |
| `Phase4_trees/gene_trees/{marker_id}.nwk` | Newick | 单棵基因树缓存（`filter` 与 `infer` 复用） |
| `Phase4_alignments/markerfinder.partition.nex` | Nexus | 分区文件 |
| `Phase5_metadata/run_config.json` | JSON | 完整参数快照与数据库版本哈希（`gtdb_markers` 或 `marker_hmm_dir`），用于复现 |
| `.markerfinder/.pipeline_state.json` | JSON | 分步运行时的显式步骤索引（分步模式使用） |
| `Phase4_intermediate/markers/{marker_id}.{faa,aln,aln.trim}` | FASTA | `--save-intermediates` 开启时保存的标记序列与比对 |
| `Phase4_intermediate/supermatrix/markerfinder.concat.fasta` | FASTA | `--save-intermediates` 开启时保存的串联比对 |
| `Phase4_intermediate/quality/checkm.tsv` | TSV | `--save-intermediates` 开启时保存的 CheckM 质量结果 |

> **阶段前缀说明**：MarkerFinder 在输出目录中实际创建 `Phase4_*`（系统发育推断与中间产物）
> 和 `Phase5_*`（报告与元数据）子目录。Phase 0～Phase 3 为内部逻辑阶段，不会作为目录输出。

---

## 资源需求

| 数据集规模 | CPU | 内存 | 预估时间 |
|-----------|-----|------|---------|
| 10 个基因组，60 个标记 | 4 | 4 GB | 15 分钟～30 分钟 |
| 50 个基因组，80 个标记 | 8 | 8 GB | 1 小时～2 小时 |
| 200 个基因组，100 个标记 | 16 | 16 GB | 3 小时～6 小时 |
| 500+ 个基因组 | 32+ | 32 GB+ | 12+ 小时 |

> 系统发育步骤是 HGT 筛查里计算密集度最高的阶段。所有候选标记均会执行该步骤。

---

## 测试

```bash
# 运行所有单元测试
pytest tests/unit -q

# 带覆盖率报告
pytest tests/ --cov=markerfinder --cov-report=html

# 运行特定测试模块
pytest tests/unit/test_marker_selection.py -v

# 环境自检（49 项检查；未获取参考数据库时为 48 项检查；只要有一项失败退出码就非 0）
markerfinder --check

# 全功能验证层：真实流水线 + 真实基因组
python validation/run_validation.py --all -n 8
```

**测试覆盖（实测值，不再是估计值）：** 测试共 88 个模块，收集 **988** 个用例，
其中 `tests/unit/` 920 个、`tests/integration/` 17 个、`tests/benchmark/` 51 个
（三个目录之和 == 徽章数字，由 `tests/unit/test_docs_numbers_are_current.py` 双向复核）。
在受支持解释器（CPython 3.10.20、3.11.15、3.12.13，均安装 ete3）下 **0 failed、0 skipped**。
用例总数以顶部徽章为准（`python3 -m pytest tests --collect-only -q`），每次提交随实跑结果刷新。
仓库内不含 CI 流水线定义（无 `.github/`），徽章仅陈述用例总数。

覆盖范围包括核心配置、模型、标记选择策略、HGT 决策引擎、GTDB-TK 标记加载、MAD 定根、
分类学解析、树与序列验证、报告生成。此外还包含三类专项测试：

- 各模块行为等价性的回归测试（`tests/unit/test_module_regressions.py`）
- must-fail 控制矩阵
- 产物级验收（`tests/integration/test_products_acceptance.py`）：在 mock 工具下真实运行流水线，
  再读回产物逐项核对

集成测试与基准测试均已有用例：`tests/integration/` 负责产物级验收，`tests/benchmark/` 提供
取数、嵌合阳性对照、度量、取证与反向消融脚本，以及对应的契约测试。仍然缺少的，是真实建树工具链
与真实基因组下的端到端数值正确性（见 MANUAL 的限制说明）。

> **Python 版本说明：** 所依赖的 ete3 当前不兼容 Python 3.13 及以上（标准库 `cgi` 已移除），
> 请使用 Python 3.10～3.12。

---

## 可复现性

MarkerFinder 每次运行后会在 `Phase5_metadata/run_config.json` 中记录完整的参数快照，包括：

- 所有配置参数
- 数据库版本哈希（SHA256）
- 软件版本
- 时间戳

复现先前运行：

```bash
markerfinder \
    -i genomes/ \
    -o output/ \
    --config previous_run/Phase5_metadata/run_config.json
```

---

## 作者

**曾子超（Zengzichao）** — [ORCID 0000-0001-6553-970X](https://orcid.org/0000-0001-6553-970X)，
[zengzichao@sjtu.edu.cn](mailto:zengzichao@sjtu.edu.cn)

---

## 引用

如果您在研究中使用了 MarkerFinder，请引用：

```bibtex
@software{markerfinder2026,
  author       = {Zengzichao},
  title        = {MarkerFinder: Adaptive HGT-Aware Phylogenomic Pipeline},
  year         = {2026},
  version      = {0.1.0},
  url          = {https://github.com/ZengZichao/MarkerFinder},
  orcid        = {0000-0001-6553-970X},
  license      = {MIT}
}
```

---

## 贡献

欢迎贡献代码。请按以下步骤操作：

1. Fork 本仓库
2. 创建功能分支 (`git checkout -b feature/amazing-feature`)
3. 提交更改 (`git commit -m 'Add amazing feature'`)
4. 推送到分支 (`git push origin feature/amazing-feature`)
5. 提交 Pull Request

请通过 [GitHub Issues](https://github.com/ZengZichao/MarkerFinder/issues) 报告 Bug 和功能需求。

---

## 联系方式

- **问题与讨论：** [GitHub Issues](https://github.com/ZengZichao/MarkerFinder/issues)
- **邮箱：** [zengzichao@sjtu.edu.cn](mailto:zengzichao@sjtu.edu.cn)
- **ORCID：** [0000-0001-6553-970X](https://orcid.org/0000-0001-6553-970X)

---

## 许可证

本项目基于 MIT 许可证发布 — 详见 [LICENSE](LICENSE) 文件。

---

## 致谢

MarkerFinder 基于以下工具和数据库构建，在此致谢：

- [MAFFT](https://mafft.cbrc.jp/) — 多序列比对软件
- [IQ-TREE3](http://www.iqtree.org/) — 高效系统发育基因组学软件
- [ASTRAL-III](https://github.com/smirarab/ASTRAL) — 最优物种树推断
- [TIGRFAM](https://www.ncbi.nlm.nih.gov/Structure/cdd/cdd.shtml) / [Pfam](https://pfam.xfam.org/) — 蛋白家族与 HMM 标记库
- [GTDB](https://gtdb.ecogenomic.org/) — 基因组分类数据库
- [CheckM](https://github.com/Ecogenomics/CheckM) — 基因组质量评估
