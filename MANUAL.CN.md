# MarkerFinder 技术手册

**自适应 HGT 感知系统发育基因组学流水线 — 技术参考手册**

版本 0.1.0

---

## 目录

1. [简介](#1-简介)
2. [安装指南](#2-安装指南)
3. [配置参考](#3-配置参考)
4. [流水线各阶段详解](#4-流水线各阶段详解)
5. [外部工具集成](#5-外部工具集成)
6. [数据库配置](#6-数据库配置)
7. [输出文件完整参考](#7-输出文件完整参考)
8. [高级用法](#8-高级用法)
9. [故障排除](#9-故障排除)
10. [算法细节](#10-算法细节)
11. [自检层产物速查](#自检层产物速查)
12. [术语表](#术语表)

---

## 1. 简介

### 1.1 什么是 MarkerFinder？

MarkerFinder 是一个系统发育基因组学流水线，自动化完成标记基因选择、水平基因转移（HGT）筛查和物种树推断。它面向原核生物基因组，为传统固定标记集方法（如 16/27/37/38 CSCG）提供了一种动态、数据驱动的替代框架，同时保留参考标记集的生物学基础。

### 1.2 设计原则

| 原则 | 说明 |
|------|------|
| **模块化** | 五大创新路径作为独立模块，构成完整流水线 |
| **自适应** | 标记集和参数阈值根据输入数据动态调整 |
| **HGT 感知** | 独立 Phase 2 基于系统发育的 HGT 筛查（Phylogenetic） |
| **双策略** | 串联法（Supermatrix）与合并法（ASTRAL）并行推断 |
| **MAG 优先** | 原生支持宏基因组组装基因组的质量感知分析 |
| **工程化** | 模块化测试套件（用例数以 README 的实测徽章为准，单一来源不在此复述）：集成层在 mock 外部工具下真跑流水线并读回产物核对验收条款，基准层含取数、嵌合阳性对照、度量、取证与反向消融脚本；两者均**不再是占位** |

### 1.3 五大创新路径

| 路径 | 模块 | 创新点 |
|------|------|--------|
| 1 | `marker_selection` | 自适应标记基因选择，4 种策略 |
| 2 | `hgt_filter` | 基于系统发育的 HGT 筛查（单系性/MAD 定根或 RF/quartet），标记分级 |
| 3 | `phylogenetic_inference` | 双策略串联+合并法，冲突检测 |
| 4 | `report_generator` | 静态 HTML 报告 + 纯文本输出 |
| 5 | `mag_optimization` | MAG 质量感知自适应参数 |
| — | `ortholog_resolver` | 多拷贝基因直系同源物解析（BBH、图聚类） |

---

## 2. 安装指南

### 2.1 系统要求

- **操作系统：** Linux（推荐）、macOS
- **Python：** >= 3.10
- **磁盘空间：** 自包含模式（系统发育 HGT 筛查 + 小体积标记 HMM，几 MB）仅需约 1 GB～12 GB，工具本体约 500 MB。HGT 筛查只基于系统发育，不需要任何外部大型蛋白数据库。
- **内存：** 最低 4 GB，> 50 个基因组的数据集推荐 16 GB

### 2.2 Conda 环境配置

```bash
# 创建并激活环境
conda create -n markerfinder "python>=3.10" -y
conda activate markerfinder

# 安装生物信息学工具（当前主流程跳过 BBH/图聚类，blast 不需要；diamond 为预留接口可选）
conda install -c bioconda -c conda-forge \
    hmmer>=3.3.2 \
    mafft>=7.520 \
    trimal>=1.4 \
    iqtree>=3.0.0 \
    astral-tree>=5.7.0 \
    fasttree>=2.1.11 \
    checkm-genome>=1.2 \
    ete3

# 如未来使用 OrthologResolver 接口，可额外安装 diamond
# conda install -c bioconda diamond>=2.1.0

# 安装 Python 依赖
pip install biopython>=1.81 pandas>=2.0 numpy>=1.24 \
    scipy>=1.10 rich>=13.0

# 安装 MarkerFinder
pip install -e ".[dev]"
```

### 2.3 验证安装

```bash
# 验证 MarkerFinder
markerfinder --version
# 预期输出: markerfinder <版本> git:<8位提交哈希> updated:<日期>
# 版本串由 setuptools-scm 从当前 checkout 生成，不要把 "0.1.0" 当作常量比对；
# 非 git 工作副本安装时无 git: 段。

# 验证外部工具
hmmsearch -h | head -1         # HMMER 3.3.2+
mafft --version                 # MAFFT 7.520+
trimal --version                # trimAl 1.4+
iqtree3 --version               # IQ-TREE3 3.0.0+
diamond version                 # DIAMOND 2.1.0+

# 运行单元测试（用例数见 README "测试"一节的实测数字）
pytest tests/unit -q

# 环境自检（49 项检查；未获取参考数据库时为 48 项检查；任一失败即非 0 退出）
markerfinder --check
```

### 2.4 目录结构

```
MarkerFinder/
├── markerfinder/              # 主包
│   ├── __init__.py            # 公共 API 导出
│   ├── __main__.py            # CLI 入口
│   ├── _version.py            # 版本与作者信息（单一常量）
│   ├── assertions.py          # 产物断言登记表与 adjudication
│   ├── banner.py              # 启动横幅显示
│   ├── config.py              # 配置 dataclass（10 个配置类）
│   ├── config_loader.py       # YAML/TOML/JSON 配置文件加载
│   ├── exceptions.py          # 自定义异常体系
│   ├── phases.py              # 流水线 Phase 登记（目录名、标签）
│   ├── pipeline.py            # 流水线协调器
│   ├── taxonomy.py            # 分类学名称解析
│   ├── validation.py          # 树/序列/交叉验证
│   ├── cli/                   # 参数解析、子命令、环境自检
│   ├── models/                # 核心数据模型
│   ├── modules/               # 五大创新路径模块
│   └── utils/                 # 工具函数
├── tests/                     # 测试套件（模块与用例数见 README，随提交变动）
│   ├── unit/                  # 单元测试（含 must-fail 控制矩阵与产物契约）
│   ├── integration/           # 集成测试（流水线状态恢复 + mock 工具下真跑并读回产物）
│   ├── data/                  # 小型手写输入（FASTA、Newick、分类学表）
│   ├── fixtures/              # 参考树与 must-fail 控制数据
│   └── benchmark/             # 外部复现与嵌合阳性对照
├── validation/                # 全功能验证：真实流水线 + 真实基因组
│   ├── cases/                 # 验收用例，按功能域分文件
│   ├── data/                  # 随仓库发布的夹具 + 来源记录（基因组需下载）
│   ├── results/               # 本次发布所依据的存档证据
│   └── scripts/               # 数据流水线：下载基因组/标记集/夹具
├── db/                        # expected_hashes.json（HMM profile 需自行获取）
├── docs/                      # 测试计划与实测报告（中英双语）
├── scripts/                   # 工具脚本（fetch_marker_db.py、smoketest.py 等）
├── pyproject.toml             # 包配置
├── environment.yml            # 含外部工具的 conda 环境
├── config.example.yaml        # 示例配置文件
├── LICENSE                    # MIT 许可证
├── README.md                  # 首页：指向中英文文档
├── README.CN.md / README.EN.md
└── MANUAL.CN.md / MANUAL.EN.md
```

### 2.5 分步运行（子命令模式）

MarkerFinder 支持把完整流程拆分为 4 个子命令逐步执行。每步的状态持久化在输出目录下
（`.markerfinder/.pipeline_state.json` 与 `Phase5_metadata/context.json`），
便于分段调试和复用中间结果。

| 子命令 | 对应 Phase | 说明 |
|--------|-----------|------|
| `scan` | Phase 0+1+1.5 | 质量预处理 + 标记扫描 + 序列提取（需 `-i`） |
| `filter` | Phase 2 | HGT 系统发育筛查 + 标记等级（缓存基因树到 `<output>/Phase4_trees/gene_trees/`） |
| `infer` | Phase 3 | 系统发育推断（串联 + 合并法，复用 filter 缓存的基因树） |
| `report` | Phase 4 | 生成 `Phase5_reports/` HTML/文本报告、`Phase4_trees/` 物种树、`Phase4_alignments/` 分区文件 + `Phase5_metadata/run_config.json` |

```bash
# 分步执行（状态在 <output>/.markerfinder/ 之间自动传递）
markerfinder scan     -i genomes/ -o output/ -t 8 --gtdb-markers-dir gtdbtk_out/align/marker_genes
markerfinder filter             -o output/ -t 8
markerfinder infer              -o output/ -t 8 --gene-tree-builder fasttree
markerfinder report            -o output/

# 也可将子命令放在参数之后（仍需在 scan 时提供 --gtdb-markers-dir）
markerfinder -i genomes/ -o output/ -t 8 --gtdb-markers-dir gtdbtk_out/align/marker_genes scan
```

**说明：**
- 不带子命令时（`markerfinder -i... -o...`）仍按原有一次聚合 `run` 执行完整流程，向后兼容。
- 每步会校验前置步骤是否已完成；例如直接运行 `filter` 而未先运行 `scan` 会报错并提示。
- 基因树缓存位于 `<output>/Phase4_trees/gene_trees/{marker_id}.nwk`，`filter` 与 `infer` 之间自动复用，避免重复建树。每个条目旁都有 `{marker_id}.nwk.input_sha256`，即该树是从哪些（基因组 id、序列）对构建出来的 SHA-256 指纹：仅当指纹与本次运行的序列一致时才复用，因此向已有输出目录重跑（`--force`）不会悄悄拿上一次运行的树当作本次的基因树。

---

## 3. 配置参考

### 3.1 PipelineConfig

顶层配置对象，由 CLI 参数或配置字典创建。

| 字段 | 类型 | 默认值 | 说明 |
|------|------|--------|------|
| `input_dir` | `str` | `""` | 输入基因组目录 |
| `output_dir` | `str` | `"./output"` | 输出目录 |
| `output_prefix` | `str` | `"markerfinder"` | 输出文件前缀 |
| `cpus` | `int` | `1` | 线程数 |
| `tmp_dir` | `Optional[str]` | `None` | 临时文件目录。默认 `None` = 自动生成系统临时目录下 `markerfinder-{uuid8}` 子目录，**每趟 run 隔离并在流程结束时自动清理**，避免多进程并行时 domtblout / concat 文件争抢；用 `--keep-tmp` 保留。显式指定路径（如 `--tmp-dir./tmp`）则沿用该路径且不自动清理。 |
| `mode` | `str` | `"standard"` | 分析模式 |
| `keep_tmp` | `bool` | `False` | 流程结束时保留临时目录（调试用）。 |
| `sequences_path` | `Optional[str]` | `None` | 参考序列文件路径 |
| `mol_type` | `Optional[str]` | `None` | 分子类型（`None`=自动检测，或 `DNA`/`RNA`/`protein`） |
| `skip_length_check` | `bool` | `False` | 跳过序列长度一致性检查 |
| `strip_annotations` | `bool` | `False` | 去除树 NHX 注释 |

### 3.2 SelectionConfig

控制 Phase 1 标记基因选择。

| 字段 | 类型 | 默认值 | 说明 |
|------|------|--------|------|
| `min_occupancy` | `float` | `0.75` | 最低单拷贝占有率阈值 |
| `max_markers` | `int` | `60` | 最大标记数 |
| `min_hmm_score` | `float` | `20.0` | 最低 HMM bit score；Phase 0 自适应参数可能覆盖该值（高质量数据集可升至 30.0，MAG 为主数据集可降至 15.0） |
| `strategy` | `SelectionStrategy` | `GREEDY` | 选择算法 |
| `marker_hmm_dir` | `str` | `""` | TIGRFAM/Pfam 单标记 HMM 目录路径（留空则由 `marker_db_source` 自动发现；仅在 `marker_mode="hmm"` 时使用） |
| `marker_mode` | `str` | `"gtdb_tk"` | 标记发现模式：`gtdb_tk`=消费 GTDB-TK ar53/bac120 每标记原始序列（需 `--gtdb-markers-dir`）；`hmm`=基于 TIGRFAM/Pfam 单标记 HMM 库扫描输入基因组（需 `--marker-hmm-dir` 或自动发现）。无 `denovo` |
| `marker_db_source` | `str` | `"auto"` | `marker_mode="hmm"` 且 `marker_hmm_dir` 为空时的 HMM 库定位策略：`auto`=在 `db/gtdb_markers/{ar53,bac120}` 下按检测到的输入域自动发现；`cog`/`gtdb`=仅查指定来源；`none`=不自动发现（需显式提供 `--marker-hmm-dir`） |
| `gtdb_markers_dir` | `str` | `""` | `marker_mode="gtdb_tk"` 时必填：GTDB-TK 提取的每标记 FASTA 目录（每个保守蛋白一个 `.faa`/`.fasta` 文件，文件名主干即标记 id；GTDB-TK 输出为**未比对/未修剪**的原始序列，软件会自动用 MAFFT 比对 + trimal 修剪后再建树） |
| `species_tree` | `str` | `""` | `--species-tree` 的别名（见 §8 示例）。提供 Newick 参考物种树时，HGT 系统发育步骤逐个比较标记基因树与参考物种树，计算 RF/quartet 不一致度。未提供时，若 `--taxonomy-table` 可用，该步骤改用 MAD 定根 + 单系比例。两者均不可用时跳过该步骤，标记记为 `UNKNOWN`（未筛查：保留，但不参与风险分档）。 |
| `quality_weighted` | `bool` | `True` | 按基因组质量加权占有率 |

### 3.3 HGTConfig

控制 Phase 2 HGT 筛查。

| 字段 | 类型 | 默认值 | 说明 |
|------|------|--------|------|
| `enable_phylogenetic` | `bool` | `True` | 启用系统发育步骤 |
| `phylogenetic_threshold` | `float` | `0.5` | 系统发育步骤风险阈值 |
| `monophyly_rank` | `str` | `"auto"` | `gtdb_tk` + `--taxonomy-table` 时单系比例的**测量阶元**。`auto`（默认）由基因树的分类学范围推导（scope 的下一阶元；tip 在任何阶元都不一致时取 `genus`）；显式命名阶元则就在该阶元测量，scope 仅被记录。取值：`auto`/`domain`/`phylum`/`class`/`order`/`family`/`genus`/`species`（见系统发育步骤章节）。 |
| `monophyly_threshold` | `float` | `0.5` | 单系比例下限：低于该值视为系统发育不一致（HGT 倾向），系统发育步骤风险 = `1 - 单系比例`。 |
| `adaptive_far_thresholds` | `bool` | `False` | **远缘自适应放宽**（**默认关闭**，须显式开启）。当输入跨多个高阶分类单元（含 ≥ 2 目，或属比 ≥ 0.70）时，系统发育风险评分趋于饱和，会把所有标记归入等级 3 剔除，导致流水线无标记可用。开启后自动放宽 `level2_max` 到 `level2_max_far`（默认 0.95），保留足量标记建树。CLI 开关为 `--hgt-adaptive-thresholds`（开启）/ `--no-hgt-adaptive-thresholds`（显式关闭）；YAML/JSON 配置文件键为 `hgt_adaptive_thresholds: true`。注意：`adaptive_far_thresholds` 仅是 HGTConfig 数据类字段名，**不能**直接作为配置文件键使用。 |
| `far_distance_genera_ratio` | `float` | `0.70` | 触发远缘放宽的属比下界（`不同属数 / 总基因组数`）。 |
| `far_distance_min_orders` | `int` | `2` | 触发远缘放宽的目数下界。 |
| `level2_max_far` | `float` | `0.95` | 远缘放宽模式下的 `level2_max` 上限。 |
| `level_thresholds` | `Dict` | `{level1_max: 0.25, level2_max: 0.60}` | 等级边界（常规模式）。远缘放宽模式下 `level2_max` 由 `level2_max_far` 覆盖。 |

### 3.4 PhylogeneticConfig

控制 Phase 3 系统发育推断。

| 字段 | 类型 | 默认值 | 说明 |
|------|------|--------|------|
| `cpus` | `int` | `1` | 建树线程数 |
| `ufboot_replicates` | `int` | `1000` | IQ-TREE3 UFBOOT 重复次数 |
| `fast_mode` | `bool` | `False` | 使用 IQ-TREE3 `-fast` 模式 |
| `use_fasttree` | `bool` | `True` | 使用 FastTree2 替代 IQ-TREE3（遗留别名，等价于 `gene_tree_builder: fasttree`） |
| `gene_tree_builder` | `str` | `"fasttree"` | 合并法基因树构建器：`"fasttree"`（默认，FastTree2 + WAG）或 `"iqtree"`（IQ-TREE3，更彻底） |
| `coalescent_mode` | `str` | `"post-filter"` | `off`、`post-filter`、`always` |
| `iqtree_timeout` | `int` | `3600` | IQ-TREE3 超时（秒） |
| `astral_timeout` | `int` | `1800` | ASTRAL 超时（秒） |

### 3.5 MAGConfig

控制 Phase 0 MAG 预处理。

| 字段 | 类型 | 默认值 | 说明 |
|------|------|--------|------|
| `checkm_results` | `str` | `None` | 预计算的 CheckM 结果文件 |
| `skip_checkm` | `bool` | `False` | 强制跳过 CheckM，使用默认质量估计 |
| `min_completeness` | `float` | `50.0` | 最低完整度 |
| `min_marker_coverage` | `float` | `0.3` | 物种最低标记覆盖率（**遗留字段**；CLI 请改用 `--min-occupancy`，废弃别名 `--min-marker-coverage` 仍可用并作用于标记选择占用率） |
| `quality_weighted` | `bool` | `True` | 按基因组质量加权分析 |
| `max_markers` | `int` | `150` | MAG 自适应模式最大标记数 |

### 3.6 OrthologConfig

控制多拷贝基因直系同源物解析。

| 字段 | 类型 | 默认值 | 说明 |
|------|------|--------|------|
| `bbh_evalue` | `float` | `1e-10` | BBH E-value 阈值 |
| `bbh_identity` | `float` | `30.0` | BBH 序列一致性阈值（%） |
| `bbh_coverage` | `float` | `0.6` | BBH 覆盖度阈值 |
| `graph_clustering` | `bool` | `True` | 启用图聚类策略 |
| `use_diamond` | `bool` | `True` | 使用 DIAMOND 替代 BLAST+ |

### 3.7 ReportConfig

控制报告生成。

| 字段 | 类型 | 默认值 | 说明 |
|------|------|--------|------|
| `output_dir` | `str` | `"./output"` | 输出目录 |
| `output_prefix` | `str` | `"markerfinder"` | 输出文件前缀 |
| `report_format` | `str` | `"html"` | 报告格式（当前仅支持 `html`） |

### 3.8 AlignerConfig

控制比对操作。

| 字段 | 类型 | 默认值 | 说明 |
|------|------|--------|------|
| `cpus` | `int` | `1` | 比对线程数 |

> AlignerConfig / OrthologConfig 等子配置的 `tmp_dir` 字段已由主配置 `PipelineConfig.tmp_dir` 统一驱动，默认均为 per-run 隔离临时目录，不必单独配置。

### 3.9 TaxonomyConfig

控制分类学名称解析。

| 字段 | 类型 | 默认值 | 说明 |
|------|------|--------|------|
| `taxonomy_table` | `Optional[str]` | `None` | 外部分类学表格路径 |
| `taxonomy_format` | `str` | `"table"` | 分类学格式：`table`、`embedded` |
| `taxonomy_source_priority` | `str` | `"table"` | 来源优先级：`table`、`embedded` |
| `taxonomy_delimiter_mode` | `str` | `"reverse"` | 解析模式：`reverse`、`greedy`、`segment` |
| `table_sep` | `Optional[str]` | `None` | 表格分隔符（自动检测为 `None`） |
| `ignore_malformed` | `bool` | `False` | 跳过畸形行而非终止（`True`=跳过，`False`=终止） |
| `taxonomy_levels` | `Optional[str]` | `None` | 以 `level:prefix` 形式（Format B 如 `kingdom:k__`，Format A 如 `kingdom:_k_`）在内置阶元之外扩展自定义阶元。该阶元会从分类学表中解析进映射，而不再被丢弃；单系筛查的阶元阶梯仍为 `domain..species` |

---

### 3.10 自检层新增开关（其余配置类）

| 字段 | 所属 | 默认值 | 说明 |
|------|------|--------|------|
| `scan_stability_min` | `ReportConfig` | `0.6` | 阈值扫描跨档入选标记集 Jaccard 下限；低于则输出“该结论记 inconclusive”的告警。由 `--scan-stability-min` 设置并真正传入扫描器。 |
| `require_evidence_coverage` | `ReportConfig` | `0.5` | HGT 证据覆盖率（`n_actually_measured / n_markers`）低于此值时在 `pipeline_summary.txt` 顶部显著告警（仅告警，不改判据）。 |
| `strict_assertions` | `ReportConfig` | `True` | 任一 `severity=fail` 输出断言触发时，进程以退出码 4 结束。 |
| `allow_assertion_failure` | `ReportConfig` | `False` | 显式把 FAIL 降级为记录性告警；绕过行为必须写进报告。 |
| `hgt_scan` | `ReportConfig` | `False` | `--hgt-scan`：产出 `{prefix}.threshold_scan.tsv`。只重评已存的风险值，**不重跑任何外部工具**。 |
| `composition_screen` | `ReportConfig` | `False` | 组成型/GC 诊断开关。由 `hgt_steps: "phylogenetic,composition"`（或 `--hgt-steps`）派生，**没有独立同名 CLI/配置文件键**；开启后产出 `Phase5_reports/{prefix}.composition.tsv`。 |
| `enable_composition_screen` | `HeterogeneityConfig` | `False` | 与 `hgt_steps: composition` **等效**：任一为真即开启组成诊断（仍无独立 CLI/YAML 键，只能程序内设置）。同对象的 `composition_warn_threshold` 经 `PipelineConfig.heterogeneity_config` 真正参与离群判定；`composition_outlier_metric` 取非 `rcv` 值时记 WARN 而不静默忽略（本轮 收尾提交，见改造日志） |
| `hgt_mode` | `HGTConfig` | `"risk"` | 标记判据模式：`risk`（现行加权风险分档）/`consistency`（ 标记级双框架交叉一致性，`grade == CONSISTENT` 为入选条件）/`hybrid`（**两者都要过**，：risk 判据与一致性判据同时满足）。默认 `risk` 永不自动切换；`consistency`/`hybrid` 启动前须通过 §5.6 结构化门禁（校验可观察工件，不能用标志位伪造），缺前置即拒绝启动并列出缺失项。过闸后在 Phase 3 真正执行判据并产出 `excluded_profile.tsv` 与 `consistency_grade` 列；`hybrid` 的逐标记联结裁决写入评估注记（进 `hgt_evaluation.tsv` 与决策卡）。 |
| `consistency_stringency` | `HGTConfig` | `1` | 一致性强度档 1..5，门槛逐级递减。`provisional`：尚未经外部基准校准，不得作为默认判据。 |
| `hgt_mode`（报告层可见） | `ReportConfig` | `risk` | 取 `consistency` 时，`pipeline_summary.txt` 会写明裁决来自跨框架一致性分档、加权 RF+quartet 风险分**仅用于排序**；若该模式下一致性判据无法运行（缺树/参照非法），则明确写出“ demotion NOT IN EFFECT”、并提示当前分档仍为 risk 制。 |
| `min_informative_sites` | `HGTConfig` | `0` | 逐标记 PIS（简约信息位点）下限，`0` = 关闭。高于 0 时，PIS 低于该值的标记即使风险再低也会标记为 `inconclusive`（不写进“干净”结论），并记入决策卡 `pis_grade`/`pis_floor`/`pis_measured`。不改变 `ev.level`，是否剔除属发布决策。阈值尚未校准。 |

**尚未接线（诚实登记，勿当作已可用）**

1. **`--hgt-mode hybrid` 已按规范实现**（定义的联结语义：“两者都要过”）：risk 判据（Phase 2 的 `excluded`，不变）与 `grade == CONSISTENT`（Phase 3 `screen`）同时满足才计入放行；一致性腿排除的标记写入 `excluded_profile.tsv`，逐标记裁决注记进 `hgt_evaluation.tsv` 与决策卡。**单趟语义**：本次运行的树由 risk 通过集构建，一致性腿裁决的是标记级结论与产物，不会假装做了第二趟重推断。`consistency` 模式**已接入运行路径**（同上产出）；缺任一可用树时记 `NOT EXECUTED`，该列保持 `NA`，不会借用 risk 结果冒充一致性结论。默认 `risk` 模式不进该分支。
2. `--cog-category-map <tsv>`已实现：`marker_id<TAB>functional_category` 两列，跳过注释/表头/单列行；文件缺失或无可用行时记 WARN，该列渲染为 `NA` 并加报告注记，不虚构类别值。
3. 组成诊断产物（`composition.tsv`）与 PIS 门槛依赖 `.pis` sidecar，只有走 `gtdb_tk` 基因树构建路径才生成；`hmm` 路径无 sidecar，故这两项在该路径上无产物。
4. **must-pass 基线的 `markers.ids` 仍是未填骨架**（需按研究对象补全 62 条核糖体标记 id）。门禁本身已实现：ids 为空时改为检验**全部**标记基因树并明确记一条日志；ids 已填但与数据集无交集时直接中止，不报“通过”。
5. 多宇宙因子评测按 E.3 排在本版本之外。
6. **任何代码都未读取的配置字段（诚实登记；由 `tests/unit/test_config_surface_ratchet.py` 强制：新增此类字段必须要么接线、要么在此登记）**：`PipelineConfig.verify_db` 与 `PipelineConfig.databases`（数据库哈希校验只在 `--check` 里做，见 run 路径不读这两个字段）；`MAGConfig.min_completeness`（完整度门槛走 Phase 0 自适应占居率路径）；`OrthologConfig.use_diamond` 与 `OrthologConfig.graph_clustering`（DIAMOND/图聚类直系同源解析未接线，Phase 1.5 明确记日志跳过）；`HeterogeneityConfig.min_allele_freq`/`consensus_threshold`/`minor_allele_threshold`/`identity_threshold`（仅供未接入主流水线的 `MAGHeterogeneityHandler` 预留，主流程每标记只取最优命中）。
   注：`HeterogeneityConfig.enable_composition_screen` 与 `composition_warn_threshold` 已不再是死字段——前者现在与 `hgt_steps: composition` 等效（任一为真即开启组成诊断），后者经 `PipelineConfig.heterogeneity_config` 真正传到 `composition.tsv` 的离群判定；`composition_outlier_metric` 取非 `rcv` 值时会记 WARN 说明只实现了 rcv，不再静默忽略。

7. **`--check --db-dir <目录>` 现在按你指定的目录做哈希校验**。此前 `--check` 内部调用哈希检查时不带参数，
   永远校验仓库默认的 `db/`，即 `--check --db-dir X` 会给出一个针对**错误对象**的结论。
   未显式传 `--db-dir` 时行为不变（仍看仓库 `db/`）。另注意 `--check` 在读取 `--config` 之前执行，
   因此想让它校验指定库请把 `--db-dir` 显式写在命令行上，不要只写在配置文件里。

---

## 4. 流水线各阶段详解

### 4.1 Phase 0：质量感知预处理

**目的：** 评估基因组质量，计算自适应参数。

**流程：**
1. **评估质量来源（`Quality source`）：**
   - **`.faa` 蛋白序列输入主动跳过 CheckM**（蛋白上预测不可靠），直接使用默认质量估计（`default_no_nucleotide`）；
   - 用户提供 `--checkm-results`（基于对应 `.fna` 预先计算）→ 预计算结果（`precomputed`）；
   - **显式 `--skip-checkm`** → 跳过 CheckM，使用默认质量估计（`skipped`）。
   - 注：当前版本强制要求 `.faa` 输入，不再支持 `.fna-only`；因此真实 CheckM 现场运行（`checkm`）与 `default_no_checkm` 来源已不存在。
2. 将基因组分层（高质量/中等/低质量/污染）
3. 根据质量分布计算自适应参数

**质量分层标准（MIMAG）：**

| 层级 | 完整度 | 污染度 |
|------|--------|--------|
| 高质量 | > 90% | < 5% |
| 中等质量 | 70%～90% | < 10% |
| 低质量 | < 70% | 任意 |
| 污染 | 任意 | >= 10% |

**自适应参数规则：**

| 高质量比例 | min_occupancy | max_markers | min_hmm_score |
|-----------|---------------|-------------|---------------|
| > 70% | 0.75 | 60 | 30.0 |
| 30%～70% | 0.55 | 80 | 20.0 |
| < 30% | 0.35 | 150 | 15.0 |

### 4.2 Phase 1：自适应标记基因选择

**目的：** 为输入数据集动态选择最优标记基因集。

**流程：**
1. **标记加载/扫描：**
   - `gtdb_tk` 模式（默认）：直接读取 `--gtdb-markers-dir` 下 GTDB-TK 提取的 ar53/bac120 每标记 FASTA，每个文件即一个标记，文件名主干为标记 id。
   - `hmm` 模式：对每个基因组运行 `hmmsearch`，扫描 `--marker-hmm-dir`（或自动发现的 `db/gtdb_markers/{ar53,bac120}`）下的 TIGRFAM/Pfam 单标记 HMM 文件。
2. **状态判定：** 对每个 (基因组, 标记) 对：
   - 0 个有效命中（score >= 阈值）→ `ABSENT`（缺失）
   - 1 个有效命中 → `SINGLE_COPY`（单拷贝）
   - 2+ 个有效命中 → `MULTI_COPY`（多拷贝）
3. **占有率计算：**
   - `occupancy(marker) = n_single_copy / n_total_genomes`
   - 质量加权变体：`Σ(completeness × [state==SINGLE]) / Σ(completeness)`
4. **标记选择：** 使用配置的策略选择最优子集

> **标记发现模式：** 默认 `marker_mode="gtdb_tk"`，直接消费 GTDB-TK 输出，**不做基因组内部相互比对**，既不会用 HMMER 重新搜索标记，也不会用 DIAMOND 做基因组之间的 BBH/自比对；标记严格取自 GTDB-TK 的单拷贝保守蛋白。GTDB-TK 的 ar53/bac120 输出是**未对齐、未修剪**的原始序列，因此每个保守蛋白在**建树前会用 MAFFT 对齐 + trimal 修剪**（这是建基因树的必要预处理，而非额外的标记发现比对）；之后逐标记建树，依据与物种树/分类学的不一致度判定 HGT，最终仅剔除高 HGT 标记、保留其余全部保守蛋白。
> `hmm` 模式用于没有 GTDB-TK 预运行的场景：MarkerFinder 用 HMMER 扫描 TIGRFAM/Pfam HMM 目录，自动识别标记并提取匹配蛋白序列。`--marker-mode` 只接受 `gtdb_tk` 与 `hmm`，无 `denovo`。

**选择策略：**

| 策略 | 算法 | 适用场景 |
|------|------|---------|
| `GREEDY` | 按单拷贝占有率降序取前 N | 标准数据集 |
| `INFO_MAX` | 带占有率重叠惩罚的贪心选择 | 多样性数据集 |
| `RATE_BALANCED` | 综合占有率 + 多拷贝比例偏离度打分 | 深层系统发育 |
| `SPARSE_OPTIMIZED` | 只要任一基因组有命中即纳入，按占有率排序 | MAG 数据集 |

**候选标记池：** `gtdb_tk` 模式下标记池等于 GTDB-TK ar53/bac120 保守蛋白集合；`hmm` 模式下标记池由 `--marker-hmm-dir` 中 TIGRFAM/Pfam HMM 文件决定。

### 4.2.5 Phase 1.5：标记序列提取与直系同源物解析

在标记基因集确定后，系统从各基因组的蛋白质 FASTA 文件中提取对应的标记基因序列。`OrthologResolver` 模块已实现基于 BBH、图聚类和长度一致性的三级直系同源物解析接口，但当前版本主流程 `_resolve_orthologs` 直接跳过该模块：在 `gtdb_tk` 与 `hmm` 模式下，主流程把输入标记当作单拷贝或已选中的最佳命中，不再进行额外的 DIAMOND BBH 解析。多拷贝情形目前采用最佳 bitscore 去重。该接口为后续版本预留，当前主流程使用最佳命中去重策略。

### 4.3 Phase 2：基于系统发育的 HGT 筛查

**目的：** 识别并按 HGT 风险对标记基因分级。

Phase 2 仅执行基于系统发育的 HGT 筛查。对每个候选标记基因，先使用 FastTree2（默认 `-wag`）构建基因树，再通过外部参照检测基因树拓扑异常。

**系统发育步骤（Phylogenetic Step）**

检测基因树与参照之间的 HGT 信号。该步骤需要外部参照才能启用，当前支持两种参照来源，
可单独提供，也可同时提供：

**(a) 参考物种树（`--species-tree`，`--tree` 为其别名）：** 提供 Newick 参考物种树时，比较每个标记基因树与参考物种树的不一致性，使用 RF/quartet 指标：

| 指标 | 权重 | 方法 |
|------|------|------|
| 归一化 RF 距离 | 50% | Robinson-Foulds 距离 / 最大可能 RF |
| Quartet 一致性 | 50% | 拓扑匹配的 quartet 比例 |

**(b) 分类学表（`--taxonomy-table`）：** 不使用物种树比较，而是对每个标记基因树先用 **MAD（Minimal Ancestor Deviation, Tria 2017）** 定根，**再统计单系比例（monophyly proportion）**。测量所用的阶元由 `--monophyly-rank` 决定：

- **`auto`（默认）**：先判断整棵树的所有 tip 在哪个阶元上完全一致（即“树的 scope”），然后在**紧邻其下的那个阶元**上统计单系性——域树查门、门树查纲、纲树查目、目树查科、科树查属、属树查种。例如某标记基因树的所有 tip 同属一个门，则 scope=门，自动在纲阶元统计单系比例。若 tip 在任何阶元都不一致（如跨多个域），`auto` 在 `genus` 上测量。
- **显式命名阶元**：就在该阶元上统计；检测到的 scope 仅被记录，不会被采用。
- **哪些分类单元会被计入**：仅当该分类单元在树上拥有 ≥ 2 个代表 tip，**且**树上不属它的 tip 也 ≥ 2 时，才计入分母。另一侧只有 1 个 tip 时，“该分类单元成单系”是无法被任何拓扑证伪的断言，因此按“不可测量”上报而不计入；分母为空则该标记跳过本筛查，并记录原因。
- **什么算单系**：基因树必须携带 `{分类单元 | 其余}` 这个 split，即该分类单元**或其补集**为 clade。因此结论与定根位置无关，不依赖 MAD 把根放在哪里。
- 单系比例 = 树上携带其 split 的分类单元数 / 可测分类单元总数；
- 单系比例低 ⇒ 基因树与分类学层级冲突 ⇒ HGT 倾向高。系统发育步骤风险 = `1 - 单系比例`，低于 `monophyly_threshold`（默认 0.5）即标记为系统发育可疑（高 HGT 风险，等级 3 剔除）。
- **阶元回退**：若正在测量的阶元下没有任何可测分类单元（例如测试集每个属仅 1 个物种），MarkerFinder 会转向其他阶元直到找到可测的那个，并在日志中记录实际使用的阶元（`Phase5_reports/{prefix}.hgt_evaluation.tsv` 的 `rank_used` 列）。当**显式命名**的阶元不得不移动时，所得比例与为该阶元设定的阈值不可比：该标记作为跨阶元比较报为 UNKNOWN，而不参与分级。仅当没有任何阶元提供可测分类单元时才真正跳过该标记的系统发育步骤。

> 注：单系性本身与定根位置无关；MAD 定根是该层统一采用的定根方法，使每个标记基因树有一个一致的、无偏的根。基因树 tip 标签使用基因组 id，因此必须与 `--taxonomy-table` 的首列一致。

**风险评分**

由于 Phase 2 仅保留系统发育步骤，HGT 风险直接等于系统发育步骤风险：

```
HGT_risk = phylogenetic
```

其中 `phylogenetic` 为系统发育步骤的综合风险（RF/quartet 不一致度或 `1 - 单系比例`）。

**等级分级（常规模式）：**

| 等级 | 风险范围 | 动作 |
|------|---------|------|
| 等级 1（清洁） | < 0.25 | 优先使用 |
| 等级 2（可疑） | 0.25～0.60 | 降权使用 |
| 等级 3（排除） | >= 0.60 | 从分析中移除 |

**远缘自适应放宽（默认关闭，须显式开启）：**

当输入跨多个高阶分类单元（含 ≥ 2 目，或属比 ≥ 0.70）时，系统发育风险评分趋于饱和，
会把所有标记归入等级 3 剔除，导致流水线无标记可用。开启后，MarkerFinder 把 `level2_max`
放宽到 `level2_max_far`（默认 0.95），保留足量标记建树。CLI 开关为
`--hgt-adaptive-thresholds`，详见 §3.3 HGTConfig。

### 4.4 Phase 3：双策略系统发育推断

**模式 A：串联法（Supermatrix）**

1. 从每个基因组提取标记蛋白序列
2. 对每个标记运行 MAFFT `--auto` 比对
3. 用 trimAl `-automated1` 修剪
4. 串联为超矩阵，缺失标记用 gap（`-`）填充
5. 写入 Nexus 分区文件
6. 用 IQ-TREE3 + ModelFinder 建树（`-m MFP -B 1000 -bnni`）；IQ-TREE3 失败时自动回退至 FastTree2

**输出：** `Phase4_trees/markerfinder.species_tree_concat.newick`、`Phase4_alignments/markerfinder.partition.nex`

**模式 B：合并法（ASTRAL）**

1. 提取并比对标记序列（同模式 A）
2. 构建单标记基因树（默认 FastTree2 + WAG，可选 IQ-TREE3，由 `--gene-tree-builder` 选择）
3. 过滤基因树（最少 4 个叶节点）
4. 运行 ASTRAL-III 推断合并物种树

**输出：** `Phase4_trees/markerfinder.species_tree_astral.newick`（仅 ASTRAL-III 成功时输出）。
ASTRAL-III 未安装或运行失败时，**不再回退到单棵基因树或 consensus**，也不输出空树文件，
只在日志中记录 ERROR。此时 `Phase4_trees/markerfinder.gene_trees.newick` 基因树集合仍然保留。
实际建树来源见 `Phase5_reports/markerfinder.pipeline_summary.txt` 的 `Species tree source` 字段。

**冲突检测**

两种模式都运行时，MarkerFinder 计算串联树与合并树之间的归一化 Robinson-Foulds 距离，
以及 quartet 一致性，据此给出拓扑冲突程度和推荐策略。当前实现**不**对单个冲突作成因分类
（HGT / ILS / 方法偏差）。

**树推荐逻辑：**

| 归一化 RF | 推荐 |
|----------|------|
| < 0.1 | 高置信 — 两种方法一致 |
| 0.1～0.3 | 中等 — 优先使用合并法（对 ILS 更稳健） |
| >= 0.3 | 低 — 强烈冲突，需检查基因树 |

### 4.5 Phase 4：报告生成

生成两类输出：
1. **静态 HTML 报告** — 自包含 HTML 摘要表，保存到 `Phase5_reports/markerfinder.report.html`（当前不支持 Plotly 交互式可视化或 PDF）
2. **纯文本输出** — TSV、Newick、Nexus 格式文件 + `Phase5_metadata/run_config.json` 参数快照（详见第 7 节）

---

## 5. 外部工具集成

### 5.1 HMMER（hmmsearch）

- **用途：** 在 TIGRFAM/Pfam HMM 配置文件中搜索蛋白序列（`hmm` 模式）
- **命令：** `hmmsearch --noali --cpu N --domtblout <out> <hmm_db> <fasta>`
- **输出格式：** DOMTABLOUT（域表）
- **解析的关键字段：** target name、query name、bit score、E-value、比对坐标（hmm_from、hmm_to、ali_from、ali_to）

### 5.2 MAFFT

- **用途：** 多序列比对
- **命令：** `mafft --auto --quiet --thread N <input>`
- **模式：** 根据输入大小自动选择算法

### 5.3 trimAl

- **用途：** 自动比对修剪
- **命令：** `trimal -automated1 -in <input> -out <output>`
- **方法：** Automated1 启发式（基于比对统计的 gap 阈值）

### 5.4 IQ-TREE3

- **用途：** 最大似然树推断
- **命令：** `iqtree3 -s <aln> -m MFP -B 1000 -bnni -nt N -pre <prefix>`
- **特性：** ModelFinder Plus 自动模型选择，UFBOOT 分支支持度

### 5.5 ASTRAL-III

- **用途：** 从基因树推断合并物种树
- **命令：** `astral -i <gene_trees> -o <output> -t 2`
- **输出：** 带局部后验概率（LPP）的物种树

### 5.6 FastTree2

- **用途：** 快速近似 ML 树（基因树）
- **命令：** `FastTree -wag -quiet -out <output> <alignment>`

### 5.7 DIAMOND

- **用途：** 快速蛋白同源搜索。用于 `OrthologResolver` 的 BBH 直系同源物解析接口（当前主流程已跳过，接口预留）
- **命令：** `diamond blastp --query <faa> --db <db> --evalue <e> --max-target-seqs N --outfmt 6...`

### 5.8 CheckM

- **用途：** 基因组质量评估（完整度、污染度）
- **命令：** `checkm lineage_wf -t N -x fna --tablename checkm.tsv <input_dir> <output_dir>`
- **输出：** 包含 Bin Id、Completeness、Contamination 列的 TSV
- **输入建议：** 当前版本强制要求 `.faa` 蛋白序列输入，`.fna-only` 不再支持。**蛋白序列输入时 MarkerFinder 主动跳过 CheckM**（蛋白上完整度预测不可靠），直接用默认质量估计（完整度 98%/污染 1%/得分 93%，MAG 来源则 75%/5%/50%）；报告 `Quality source = default_no_nucleotide`。如需使用预计算的 CheckM 质量，可通过 `--checkm-results` 提供基于对应 `.fna` 的结果文件，此时 `Quality source = precomputed`。
  - 使用 `--skip-checkm` 可显式跳过 CheckM，此时 `Quality source = skipped`。
  - 注：由于输入已强制为 `.faa`，真实 CheckM 现场运行（`checkm`）与 `default_no_checkm` 来源已不再出现。

---

## 6. 数据库配置

### 6.1 `hmm` 模式：TIGRFAM/Pfam 单标记 HMM 目录

`hmm` 模式下，MarkerFinder 需要一个包含**单个标记一个 HMM 文件**的目录（如 GTDB-Tk-214-Markers 的 `pfam/individual_hmms/` 和 `tigrfam/individual_hmms/`）。文件名主干即标记 id，后缀为 `.HMM` 或 `.hmm`。代码会自动递归搜索并排除合并库文件（如 `Pfam-A.hmm`）。

```bash
# 方式一：显式提供 per-marker HMM 目录
markerfinder -i genomes/ -o output/ --marker-mode hmm \
  --marker-hmm-dir GTDB-Tk-214-Markers/pfam/individual_hmms

# 方式二：拼装 db/gtdb_markers/{ar53,bac120} 供 marker_db_source=auto 自动发现。
# 这些 profile 属第三方模型，本仓库不入库，请从本地的 GTDB-Tk 安装（或你已持有的文件）拼装：
python scripts/fetch_marker_db.py --auto
python scripts/fetch_marker_db.py --set bac120 --from-hmm .../gtdbtk_bac120.a.hmm
python scripts/fetch_marker_db.py --set ar53  --from-dir .../ar53_marker_genes
# 加 --record 会顺便把拼装后的目录哈希写入 db/expected_hashes.json，
# 那正是 `markerfinder --check` 比对的目标。

# 验证
hmmstat db/gtdb_markers/bac120/*.HMM | head
markerfinder --check                       # 报告计算值与期望值
```

### 6.2 `gtdb_tk` 模式：GTDB-TK 每标记 FASTA（必需）

该模式不读 `db/`：它消费的是 `gtdbtk align` 针对**你的输入基因组**输出的每标记 FASTA，
由 `--gtdb-markers-dir` 传入。

```bash
# 先对输入基因组跑一次 GTDB-Tk（对齐与提取都由它完成）
gtdbtk align --genome_dir genomes/ --out_dir gtdbtk_out/ --cpus 8

# MarkerFinder 读 gtdbtk_out/align/marker_genes/{bac120,ar53}/<marker>.faa
markerfinder -i genomes/ -o output/ --marker-mode gtdb_tk \
  --gtdb-markers-dir gtdbtk_out/align/marker_genes
```

每个标记一个 FASTA，文件名主干即标记 id；对齐所用的 GTDB-Tk 版本会逐次记录（见 6.4）。

### 6.3 gtdb_tk 模式示例（系统发育步骤走 MAD 定根 + 单系比例）

```bash
# gtdb_tk 模式 + 分类学表（系统发育步骤走 MAD 定根 + 单系比例，无需物种树）
markerfinder -i genomes/ -o output/ \
  --marker-mode gtdb_tk \
  --gtdb-markers-dir gtdb_out/ \
  --taxonomy-table gtdb.summary.tsv \
  --monophyly-rank auto \
  --monophyly-threshold 0.5
```

### 6.4 版本锁定

每次运行后，MarkerFinder 在 `Phase5_metadata/run_config.json` 中写入完整参数快照和数据库版本哈希，确保分析可复现：

路径字段记录的是**本次运行实际使用的绝对路径**。相对输出目录的值在事后无法解释（它是相对于谁而言的？），
而若按当前工作目录去解，重放就会写到没人要求的地方。因此用
`--config <output>/Phase5_metadata/run_config.json` 回喂时，会恢复 marker 目录、taxonomy 表与风险带，
会恢复用户**请求的**标记预算与占居率下限（而不是 Phase 0 自适应后的值），并丢弃上一次运行的临时目录
（暂存空间属于正在重放的那次运行）。命令行无法表达的字段（如 `database_versions`、`databases`）
会以警告列出，而不是静默忽略。

```json
{
  "markerfinder_version": "0.1.0",
  "timestamp": "2025-07-02T10:30:00",
  "run_duration_seconds": 1234.5,
  "parameters": { "mode": "standard", "min_occupancy": 0.75, ... },
  "database_versions": {
    "gtdb_markers": {"name": "gtdb_markers", "path": "db/gtdb_markers/bac120", "hash": "abc123...", "exists": true}
  }
}
```

复现先前运行：

```bash
markerfinder -i genomes/ -o output/ --config previous_run/Phase5_metadata/run_config.json
```

---

## 7. 输出文件完整参考

输出目录按运行阶段分为以下子目录（名称带 `PhaseN_` 前缀）：

- `Phase5_reports/` — HTML/文本报告与 TSV 汇总表
- `Phase4_trees/` — 物种树与基因树集合
- `Phase4_trees/gene_trees/` — 单棵基因树缓存（每个标记一个 `{marker_id}.nwk`，并附其溯源指纹 `{marker_id}.nwk.input_sha256`）
- `Phase4_alignments/` — Nexus 分区文件
- `Phase5_metadata/` — `run_config.json`、`context.json`
- `.markerfinder/` — 分步运行内部状态 `.pipeline_state.json`

### 7.1 HTML 报告（`Phase5_reports/markerfinder.report.html`）

自包含静态 HTML 文件，包含：

- 执行摘要（基因组数、标记数、运行时间）
- HGT 等级分布（L1/L2/L3 计数）
- 各输出文件路径说明

> 当前仅支持 HTML 静态报告，不支持 Plotly 交互式可视化或 PDF。

### 7.2 标记汇总（`Phase5_reports/markerfinder.marker_summary.tsv`）

```tsv
marker_id	occupancy_score	marker_quality_score	marker_quality_level
COG0049	0.9500	0.8234	level_1
COG0085	0.8800	0.7856	level_1
```

### 7.3 HGT 评估（`Phase5_reports/markerfinder.hgt_evaluation.tsv`）

```tsv
marker_id	overall_risk	hgt_risk_level	hgt_evidence_confidence
COG0049	0.1200	level_1	high
```

### 7.4 物种树（`Phase4_trees/markerfinder.species_tree_concat.newick`）

串联法物种树，标准 Newick 格式，带分支支持度值。

### 7.5 分区文件（`Phase4_alignments/markerfinder.partition.nex`）

```nexus
#nexus
begin sets;
  charset COG0049 = 1-300;
  charset COG0085 = 301-600;
end;
```

### 7.6 流水线摘要（`Phase5_reports/markerfinder.pipeline_summary.txt`）

人类可读的文本文件，包含:

> **未测量声明**：本次运行中因外部工具或依赖不可用而**没有测到**的指标，会在摘要顶部逐条列出（指标名 + 原因 + 受影响标记），同样的内容也以 WARNING 级写入日志。这些指标在表格里一律渲染为 `NA`，绝不以 0 或 0.5 这类中性值出现；覆盖率数字（`Evidence coverage`）只回答“有多少标记真正测过”，而这份声明回答“缺的是哪一项、为什么缺”——两者不可互相替代。


```
============================================================
MarkerFinder Pipeline Summary
============================================================

Genomes: 6
Markers: 60
HGT L1: 0 L2: 30 L3: 30
Species tree source: concat
Quality source: default_no_nucleotide
Runtime: 57.8s
```

- `Species tree source`: 取值 `concat`（串联法）/ `astral`（合并法 ASTRAL-III 成功）/ `none`（ASTRAL-III 失败或未跑 coalescent）。当前实现不会在 ASTRAL-III 失败时生成 `consensus` 或 `first_gene_tree` 文件。
- `Quality source`: 取值 `precomputed`(用户提供.fna 质量) / `default_no_nucleotide`(蛋白.faa 输入主动回退) / `skipped`(`--skip-checkm` 显式跳过)。当前版本强制 `.faa` 输入，真实 CheckM 现场运行来源 `checkm` 与 `default_no_checkm` 已不再出现。

### 7.7 合并法物种树（`Phase4_trees/markerfinder.species_tree_{source}.newick`）

- **串联法物种树** `Phase4_trees/markerfinder.species_tree_concat.newick`：Phase 3 模式 A（Supermatrix）输出，默认以 IQ-TREE3 推断，失败时自动回退至 FastTree2，**默认始终生成**。
- **合并法物种树** `Phase4_trees/markerfinder.species_tree_astral.newick`：仅在启用 `--coalescent-mode` 且 ASTRAL-III 正常退出时输出；失败时不输出空树文件，仅记录 ERROR。
- `Quality source` 与 `Species tree source` 字段同时写入 `Phase5_reports/markerfinder.pipeline_summary.txt`（见上）。

### 7.8 基因树集合（`Phase4_trees/markerfinder.gene_trees.newick`）

各标记基因的 Newick 树集合，标准 multi-newick 格式（每行一棵树），无 `>gene_id` 分隔头。仅在启用合并法推断时生成。单棵基因树缓存位于 `Phase4_trees/gene_trees/{marker_id}.nwk`。

### 7.9 运行配置快照（`Phase5_metadata/run_config.json`）

JSON 格式的完整参数快照，包含所有配置参数、数据库版本 SHA256 哈希（主要校验 `gtdb_markers` 或 `marker_hmm_dir`）、软件版本和时间戳，用于复现分析（详见 6.4 节）。

---

## 8. 高级用法

### 8.1 使用预计算的 CheckM 结果 / 跳过 CheckM

```bash
# 预计算 CheckM 结果
markerfinder -i genomes/ -o output/ -t 8 \
    --checkm-results /path/to/quality_report.tsv

# 显式跳过 CheckM（例如快速测试）
markerfinder -i genomes/ -o output/ -t 8 \
    --skip-checkm
```

CheckM 结果文件应为 TSV 格式，包含列：`Bin Id`、`Completeness`、`Contamination`。当前版本强制要求输入目录包含 `.faa` 文件，不再支持 `.fna-only`；`--checkm-results` 仅用于导入基于历史 `.fna` 的预计算质量。
使用 `--skip-checkm` 时，`pipeline_summary.txt` 中 `Quality source` 会显示为 `skipped`。

### 8.2 自定义 HGT 阈值与远缘放宽

```bash
# 常规模式：level1_max 取 --hgt-threshold，常规 level2_max = 0.60
markerfinder -i genomes/ -o output/ -t 8 \
    --hgt-threshold 0.30

# 远缘数据集显式开启自适应放宽（把 level2_max 上限放宽到 0.95；默认关闭）
markerfinder -i genomes/ -o output/ -t 8 \
    --hgt-adaptive-thresholds
```

> 远缘放宽**默认关闭**：需要时请显式传 `--hgt-adaptive-thresholds`（或在配置文件写 `hgt_adaptive_thresholds: true`），MarkerFinder 才会对远缘数据集（含 ≥ 2 目，或属比 ≥ 0.70）把 `level2_max` 放宽到 0.95，缓解系统发育风险饱和致使 60/60 全部等级 3 剔除的问题（详见 §3.3 HGTConfig 表后方的远缘放宽一节）。

### 8.3 大规模数据集优化

```bash
markerfinder -i genomes/ -o output/ -t 32 \
    --mode standard \
    --gene-tree-builder fasttree \
    --coalescent-mode post-filter
```

关键优化：

- `--gene-tree-builder fasttree`：使用 FastTree2 + WAG 构建基因树（快 10 倍，默认）
- `--coalescent-mode post-filter`：仅在 HGT 过滤后运行合并法

### 8.4 跳过合并法推断

```bash
markerfinder -i genomes/ -o output/ -t 8 --coalescent-mode off
```

### 8.5 使用配置文件

通过 `--config` 参数加载 YAML、TOML 或 JSON 配置文件，减少重复的命令行输入：

```bash
markerfinder -i genomes/ -o output/ --config config.yaml
```

配置文件示例（`config.yaml`）：

```yaml
threads: 8
mode: standard
hgt_steps: "phylogenetic"
hgt_threshold: 0.25
coalescent_mode: post-filter
ufboot: 1000
```

**优先级规则：** CLI 参数 > 配置文件 > 内置默认值。即命令行中显式指定的参数始终覆盖配置文件中的值。

### 8.6 输出文件管理

```bash
# 强制覆盖已有输出
markerfinder -i genomes/ -o output/ --force

# 跳过已有文件（不覆盖）
markerfinder -i genomes/ -o output/ --no-clobber
```

默认行为：如果输出目录已包含文件，MarkerFinder 将报错退出，提示使用 `--force` 或 `--no-clobber`。

### 8.7 临时目录管理

MarkerFinder 默认在系统临时目录下自动生成 `markerfinder-{uuid8}` 作为临时根目录，每趟 run 隔离，并在流程结束时**自动清理**。

```bash
# 保留临时目录用于调试( FAILED 时可查看 .domtblout/.faa/.aln/ )  
markerfinder -i genomes/ -o output/ --keep-tmp

# 显式指定临时目录( 不自动清理 )
markerfinder -i genomes/ -o output/ --tmp-dir ./my_tmp
```

### 8.8 与原版 MarkerFinder 的向后兼容

从原版 MarkerFinder 脚本迁移的用户：

- 使用 `--mode conservative` 获取最接近原版固定标记集的行为
- 原版 16、27、37 和 38 号 CSCG 集已包含在候选 COG 池中

---

## 9. 故障排除

### 9.1 常见错误

| 错误 | 原因 | 解决方案 |
|------|------|---------|
| `NoMarkerAvailableError` | 没有标记满足占有率阈值 | 使用 `--mode mag_adaptive` 或 `--mode expanded` |
| `hmmsearch not found` | HMMER 未安装 | `conda install -c bioconda hmmer` |
| `HMM database not found` | 数据库路径未配置 | 设置 `--marker-hmm-dir` 或将 TIGRFAM/Pfam HMM 文件放入 `db/gtdb_markers/{ar53,bac120}`；或改用 `gtdb_tk` 模式（`--gtdb-markers-dir`） |
| `gtdb-markers-dir not provided` | `gtdb_tk` 模式缺少必填参数 | 提供 GTDB-TK 提取的每标记 FASTA 目录 |
| `IQ-TREE3 failed` | 内存不足或比对问题 | 使用 `--gene-tree-builder fasttree` 或减少 `--ufboot` |
| `ASTRAL failed` | 基因树 < 4 个物种 | 确保输入至少 4 个基因组 |
| `No genomes found` | 目录错误或文件扩展名不对 | 确保输入目录中存在 `.faa`/`.fasta`/`.fa` 蛋白序列文件；当前版本不再支持 `.fna-only` |

### 9.2 性能问题

| 症状 | 诊断 | 解决方案 |
|------|------|---------|
| Phase 1 慢 | hmmsearch 瓶颈 | 增加 `--threads`，SSD 用于 tmp |
| Phase 2 系统发育步骤慢 | 每个标记运行 IQ-TREE3 | 使用 `--gene-tree-builder fasttree`（默认）构建基因树，或减少同时分析的标记数 |
| Phase 3 慢 | 大超矩阵 | 使用 `--gene-tree-builder fasttree` |
| 高内存 | 大比对 | 减少 `--ufboot`，使用 `--gene-tree-builder fasttree` |

### 9.3 降级模式行为

当数据库不可用时，MarkerFinder 优雅降级：

| 缺失资源 | Phase 1 | Phase 2 | Phase 3 |
|---------|---------|---------|---------|
| HMM 数据库 | 全 ABSENT → 报错（可改用 `gtdb_tk` 模式或显式提供 `--marker-hmm-dir`） | 系统发育步骤跳过（无参照） | 空比对 |
| CheckM | 默认估计值 | — | — |

### 9.4 退出码

MarkerFinder 使用标准化退出码，便于脚本和工作流集成：

| 退出码 | 含义 | 典型场景 |
|--------|------|---------|
| `0` | 成功 | 流水线正常完成 |
| `1` | 运行时错误 | 外部工具失败、内存不足、数据格式错误 |
| `2` | 参数错误 | 无效的命令行参数、缺少必填参数 |
| `3` | 数据错误 | 输入文件不存在、格式无法识别、树验证失败、分类表含不可解析行 |
| `4` | **输出断言失败 / must-pass 未过** | 任一 `severity=fail` 的输出数值断言触发，或分类学 must-pass 对照集有“必须成立”关系不成立。用 `--allow-assertion-failure` 可显式降级为告警，绕过行为会写进报告 |
| `5` | **判定不确定（inconclusive）** | 流程跑完但推荐层交不出确定树：冲突度量不可用时 `recommend_tree` 返回 `recommended_tree=None` / `confidence="inconclusive"`，不再“照常交出一棵树” |
| `130` | 用户中断 | 用户按 Ctrl+C 终止（SIGINT） |

---

## 10. 算法细节

### 10.1 Robinson-Foulds 距离

```
RF(T1, T2) = 不同二分分割的数量
normalized_RF = RF / (2 × (n - 3))
```

- 0.0：拓扑完全相同
- 0.3～0.5：中等不一致
- > 0.5：强烈不一致

### 10.2 标记质量评分

```
score = 0.25 × hmm_norm + 0.25 × occupancy + 0.10 × length_norm 
      + 0.20 × informativeness + 0.20 × (1 - hgt_risk)
```

| 分数范围 | 等级 |
|---------|------|
| >= 0.8 | 等级 1（高质量） |
| 0.5～0.8 | 等级 2（中等质量） |
| < 0.5 | 等级 3（低质量） |

---

## 自检层产物速查

| 产物 | 内容 | 工作包 |
|------|------|--------|
| `Phase5_reports/{prefix}.mustpass.tsv` | 分类学 must-pass 门禁的 `summary`/`violation`/`not_checked` 三类行（`--taxonomy-mustpass` 开启时产出） | |
| `Phase5_reports/{prefix}.assertions.tsv` | 每条断言的 `assertion_id/name/severity/result/detail/provisional` | |
| `Phase5_reports/{prefix}.threshold_scan.tsv` | 每档 `(L1,L2,L3,UNKNOWN)` 四元组 + `n_at_boundary` + 跨档 Jaccard 矩阵 | |
| `Phase5_evidence/decision_<n>.json` | 逐标记决策卡（机器可读，`schema_version: 2`） | |
| `Phase5_reports/{prefix}.excluded_profile.tsv` | 筛除标记的 `marker_id`、`grade`、`stringency`、`reasons`、`pis`、`functional_category`、`category_note`，以及两条腿各自的 `concat_strength`/`concat_quartets` 与 `coalescent_strength`/`coalescent_quartets`（strength= 的一致性强度，quartets=该一致性由多少个四尖端组合算出；不可测的一律 `NA`，绝不用 0 顶替） | |
| `hgt_evaluation.tsv` 右侧追加列 | `detector/rank_used/n_total/n_mono/rf/rf_state/quartet/quartet_state/…/assertion_ids_fired`（旧列列序不变，`None` 渲染为 `NA`） | |
| `marker_summary.tsv` 右侧追加列 | `pis` / `effective_columns` / `consistency_grade`（列序不变，只右追加；`consistency_grade` 在未跑判据时为 `NA`）。 另列的 `hgt_risk_score`、`info_rank` 两列**尚未实现**（风险分在 `hgt_evaluation.tsv` 里，不重复入表） | |
| `pipeline_summary.txt` 顶部 | 断言 PASS/WARN/FAIL 计数、证据覆盖率、** 未测量指标声明**、far 模式生效声明 | |

退出码新增：**4** = `EXIT_ASSERTION_FAILED`（输出数值不合理，或分类学 must-pass 未通过）；**5** = `EXIT_INCONCLUSIVE`（流程跑完但推荐层交不出确定树）。`0/1/2/3/130` 含义不变。

## 术语表

本表与 的附录 A 同源。

| 术语 | 本项目用法 | 注意事项 |
|---|---|---|
| **marker / gene** | marker = 挑选出来的直系同源标记；gene = 某基因组中该 marker 的拷贝 | 新增字段一律用 `marker_id` |
| **risk（风险分）** | `overall_risk ∈ [0,1]`，由 RF 与 quartet 两信号加权合成 | 不是概率，不可解读为 “HGT 概率” |
| **incongruence（拓扑不一致）** | 基因树与参照（物种树 / 分类框架）之间的分歧本身，由 `normalized_rf` 与 `quartet_agreement` 度量；是一个**观测量** | 与 `risk` 严格区分：risk 是把不一致程度加权合成后用于**排序**的分，不是概率，也不是 HGT 结论；不一致本身还可能来自 ILS、组装碎片化或长枝吸引，因此“有不一致”不等于“发生了 HGT” |
| **consistent / concordant** | `consistent` 专指标记级双框架同侧判据；quartet agreement 指数据集级拓扑一致比例 | 不可互换，二者是不同层级的量 |
| **不可测 (NOT_MEASURABLE)** | 依赖/工具/数据缺失导致根本无法测量 | 与“不适用”区分，用户行动不同 |
| **不适用 (NOT_APPLICABLE)** | 该分量在此判定路径上本就不该存在（如单系路径不测 RF） | 不是错误，但也不该填成 0 |
| **被拒绝 (REJECTED)** | 测了，但参照对象合法性未通过 | 错误 A 的防线（错误 A 指撤稿论文 Steenwyk & King 2025, *Science* 390:751-756, doi:10.1126/science.adw9456 的合并参照缺陷；详见 README「概念来源」节） |
| **UNKNOWN** | `MarkerLevel.UNKNOWN`，HGT 层“未筛查”，标记被保留但不参与风险叙事 | 与 `inconclusive`（推荐层/一致性层）是两回事 |
| **inconclusive** | 证据存在但强度不足 ⇒ 不下结论（推荐层/一致性层枚举） | `MarkerLevel` 不新增此值 |
| **证据覆盖率** | `n_actually_measured / n_markers` | 低于 `--require-evidence-coverage` 时报告顶部警示 |
| **must-fail 控制** | 证明某检查“真的会红”的负面 fixture | 没有它，检查的存在性等价于不存在：本套件的每个负面扫描都自带一个植入的正对照 |
| **参照树** | 被用来给基因树打分的树（物种树/串联树/ASTRAL 树/MAD 定根后的自身） | 参照对象本身必须通过合法性校验 |
| **far 模式** | `adaptive_far_thresholds`：远缘数据集把 level2_max 由 0.60 放宽至 0.95 | 生效时报告顶部与受影响标记卡均留痕 |

---

*手册结束*
