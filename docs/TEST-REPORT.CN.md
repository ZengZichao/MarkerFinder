# MarkerFinder 测试报告 — 发布验收的实测结果

测试方案：[TESTING.CN.md](TESTING.CN.md)（本报告执行的正是该方案）·
适用版本：MarkerFinder 0.1.0 · [English version](TEST-REPORT.EN.md) — 英文版为主版本。

> 本报告中每一个数字都来自随包发布的一份产物，且
> `tests/unit/test_validation_report_is_current.py` 会在每次 `pytest tests` 时从这
> 些产物重新推导主指标。数字与证据不一致会导致测试失败，而不只是一句过时的话。

---

## 1. 本次运行的标识

| 字段 | 值 |
|---|---|
| 被测软件 | MarkerFinder `0.1.0` |
| 解释器 | CPython 3.11.9，conda 环境 `markerfinder` |
| 主机 / CPU | Linux x86-64，96 核（取自产物中的 `platform.platform`） |
| 验收运行 | `python validation/run_validation.py --all -n 16`（16 个 pytest worker） |
| 用例数 | 收集 153 / 执行 153 |
| 结果 | **152 passed · 0 failed · 0 errors · 1 xfailed**（预期失败，见 §5.7） |
| 墙钟耗时 | 见 `validation/results/junit.xml` 中套件节点的 `time` 属性 |
| 外部工具 | hmmsearch 3.4、mafft v7.525、trimal v1.5.rev1、FastTree 2.2.0、IQ-TREE 3.1.3、ASTRAL 5.7.8、CheckM 1.2.2、ncbi-datasets 18.38.0 |
| 归档产物 | `validation/results/{junit.xml,metrics.json,capability_matrix.tsv,environment.txt,chimera_metrics.json}` |

完整环境清单——解释器路径、每个声明依赖的实际版本、每个外部工具自报的版本与路径、
本次运行使用的 `TMPDIR`——见 `validation/results/environment.txt`；它由
`validation/scripts/04_record_environment.py` 通过**执行命令**生成，而不是抄自文字。

## 2. 两层测试

| 层 | 用例数 | 结果 | 能证明什么 |
|---|---|---|---|
| `tests/`（单元 + 集成 + 基准驱动，外部工具 mock） | 收集 988 | 988 passed · 0 skipped · 0 failed | 算术、schema、接线、拒绝路径 |
| `validation/`（真实流水线、真实基因组、真实工具） | 153 | 152 passed · 1 xfailed · 0 failed | 结论：筛查能否区分、门禁能否触发、植入的转移能否检出 |

快速层不需要外部二进制、也不下载数据；验收层两者都需要，它是发布的验收证据。
`validation/results/capability_matrix.tsv` 记录了哪个用例覆盖哪个被声称的能力。

## 3. 功能面覆盖

功能面由代码枚举，而不是靠手工维护清单（`validation/capabilities.py`）：

| 权威来源 | 条目数 | 至少被一个用例覆盖 |
|---|---|---|
| `_build_parser` 中的命令行选项（含废弃别名及其记录属性） | 68 | 68 |
| 分步子命令（`scan` / `filter` / `infer` / `report`） | 5 | 5 |
| `markerfinder.cli.constants` 中互不相同的退出码 | 7 | 7 |
| 输出产物（报告、证据、树、比对、状态） | 19 | 19 |
| 流程性质（选择、分级、溯源、确定性、失败必须大声…） | 9 | 9 |
| **能力 id 合计** | **108** | **108 — `capabilities_uncovered: []`** |

这项检查是双向的：有能力面无用例会失败；用例声明了 parser 已不再暴露的能力 id 也会
失败。`test_v90` 会在 `metrics.json` 中报出
`capabilities_required = 108, capabilities_covered = 108`；而部分收集（`-k …`）会让
该守卫**显式 skip**，而不是拿着削弱的快照报通过。

## 4. 测试数据：溯源与完整性

| 项目 | 值 | 证据 |
|---|---|---|
| 蛋白质组 | 12 个，来自 NCBI RefSeq（`datasets download genome accession … --include protein`） | `validation/data/genomes/`、`PROVENANCE.tsv` |
| lineage / 物种名 / 组装统计 / 完整度与污染率 | 全部取自 NCBI 自己的报告（`datasets summary genome accession`、`datasets summary taxonomy taxon`） | `validation/data/PROVENANCE.tsv` 的 24 列 × 12 行 |
| marker FASTA | 用 hmmsearch 扫描 `db/gtdb_markers/bac120`（由 `scripts/fetch_marker_db.py` 拼装）的 profile 后，对 GTDB-Tk 布局的**重构**；120 个被扫描 profile 的占居率全部记录 | `validation/data/markers/*/OCCUPANCY.tsv`、`OCCUPANCY_ALL.tsv` |
| 跨目阳性对照 | 植入 3 次跨目转移（受体 *Bacillus anthracis* Ames ← 供体 *Lactococcus lactis* Il1403），从被记录的候选池以 `random.Random(20260929)` 抽取 | `validation/data/markers/chimera_ground_truth_bench12.json` |
| 功能类别映射 | 取 profile 自身的 `DESC` 注释行（8/8 个 marker 均有注释）——**不是** COG 生物学 | `validation/data/cog/marker_category_map.tsv` |
| 参照树 | 由记录的 lineage 构造；所有内部节点均已解析（未解析的多分叉会让参照不可测量），并在构造时校验 | `validation/data/trees/`、`03_build_fixtures.py::_validate_reference` |
| 完整性 | 对每个随包数据文件计算 SHA-256 | `validation/data/MANIFEST.sha256`（241 个文件），由 `test_shipped_data_matches_the_manifest` 校验 |
| 物种身份 | 每个组装记录的属名，必须是该蛋白质组 FASTA 头中 `[organism]` 的多数 | `test_shipped_genomes_are_the_organisms_the_provenance_claims` |

构建该数据包时发现并修掉了两个数据缺陷：taxonomy 表原本为手工填写，把
`GCF_000008105.1` 标成 *Bacillus*（NCBI 记录为 *Salmonella enterica*——这个错误会让
MAD 定根筛查去检验一个数据里并不存在的关系）；生成的参照 Newick 曾按分类阶元逐层重复
支长（`GCF_000009045.1:0.01:0.01:0.01`），ete3 会把它读成一个名为 `:0.01` 的节点。
这两类问题现在都有构造期的拒绝式校验。

## 5. 实测到了什么

### 5.1 标记选择的旋钮确实决定结果

| 旋钮 | 差示测量 | 产物字段 |
|---|---|---|
| `--max-markers 2 / 4 / 8` | “Markers selected” = 2 / 4 / 8 | `v07_budget.selected_per_ceiling` |
| `--min-occupancy 0.50 / 0.95` | 保留 20 个 → 保留 17 个 | `v07_occupancy` |
| `--marker-preset …` | 每个具名预设都覆盖显式旋钮；`none` 不覆盖 | `v07_preset_*` |
| `--mode conservative/standard/expanded/mag_adaptive` | 都能跑完并在摘要中标出自己的模式 | `v07_mode_*` |
| `--min-marker-coverage`（废弃别名） | 与 `--min-occupancy` 选出同样 8 个 marker，并声明自己是别名 | `v07_alias` |

`--max-markers` 过去是「读取了、记录了、然后被 Phase 0 自适应预算覆盖」的死旋钮。
现在它在自适应之后生效（显式 `--marker-preset` 同理），第一行测的正是这件事。

### 5.2 系统发育 HGT 筛查能区分一致与冲突

在 8 基因组集合上以 `--monophyly-rank order` 测量（该阶元下每个 marker 有 2 个可测分类
单元）：

| 量 | 值 | 出处 |
|---|---|---|
| order 阶元下被测试的分类单元数（`n_total`） | 2 | `hgt_evaluation.tsv` |
| 单系比例分布 | 4 个被评分 marker 中出现 0.0 与 0.5 | `hgt_evaluation.tsv` |
| `--monophyly-threshold 0.0` vs `1.0` | 3 个被标记 vs 3 个被标记 | `v04_monophyly` |
| `--hgt-threshold-l2l3 0.9` vs `0.1` | 剔除 0 个 → 剔除 2 个 | `v06_l2l3` |
| `--hgt-threshold`（L1/L2 的废弃别名） | Level-1 数 4 → 2，Level-3 数不变 | `v06_alias` |
| `--hgt-adaptive-thresholds` | `level2_max` 0.6 → 0.95，并报出远缘判定；默认运行不提 | `v06_adaptive` |
| 仅排名降级（一致性判据） | `--cog-category-map` 会被写进被剔除标记的 `functional_category`；不给图时该列为 `NA` 并附说明 | `v06_cog` |

这一层的三个行为是错的并已修复，因为是验收数据把它们暴露出来的：单系判定依赖 MAD 把
根放在哪里（手册明言单系性与定根无关，而代码只检验「有根树的 clade 成员」）；补集只有
1 个 tip 的分类单元被当作可测量（7 对 1 的取样无法证伪任何断言，代码却因此给完全一致
的树打出 risk 1.0）；`--monophyly-rank` 这个文档承诺的阶元会被一条按 scope 推导的自动
规则覆盖，用户无法关闭。现在默认是 `auto`（由树推导），显式命名阶元则以之为准。

### 5.3 门禁、断言与退出码

| 机制 | 实测结果 | 退出码 |
|---|---|---|
| 健康运行下的 14 条输出断言 | 14 PASS / 0 WARN / 0 FAIL | 0 |
| `--taxonomy-mustpass` 满足 | `groups_tested=8; markers_checked=8`——门禁确实测了东西 | 0 |
| `--taxonomy-mustpass` 被违反（genus 拆分孪生表） | 运行中止，并点名违反项 | 4（`EXIT_ASSERTION_FAILED`） |
| `--taxonomy-mustpass` 指向不存在的基线 | 拒绝而不是通过 | 4 |
| `--allow-assertion-failure` 绕过 | 绕过被记录进产物 | 4 → 运行完成但绕过被标出 |
| 推荐结论不可判定 | 文档所述的退出码 5 可达且有说明 | 5 |
| 参数错误（`-t 0`、`--min-occupancy 1.5`、`-o` 指向已存在文件） | 带说明拒绝，无 traceback | 2 |

must-pass 门禁过去会把基因树缓存文件的**路径**当作 Newick 读取，于是它一个关系都没测
却打印「所有必需关系成立」。上面这条基于产物的断言就是防止该形状的错误回归的手段。

### 5.4 推断与产物

* `--gene-tree-builder fasttree` 与 `iqtree` 都能完成；support 值只在能产出它的构建器
  下出现，并且标签写清用了哪种（`v09_support`、`v09_builder`）。
* `--coalescent-mode off / post-filter / always` → ASTRAL 树分别为不存在、存在、存在
  （`v09_coalescent_*`），推荐文字跟着真实存在的树走。
* 18+1 类产物全部被打开并解析：TSV 表头、决策卡 JSON schema、Nexus 分区文件
  （`charsets` 数等于筛查**保留**的标记数，而不是 Phase 1 的选中数）、HTML 报告、
  run_config、分步状态（`v14_products`、`v14_partition`、`v14_cards`、`v13_*`）。
* 读取物种树 tip 集合时会先剥掉 FastTree 的 support 标签，因此 `72` 这类标签不会被当成
  一个分类单元（`v14_trees`）。
* 输出目录之外不写任何文件；并发运行之间的暂存目录互相隔离
  （`v13_tmp`、`v13_tmp_isolation`、`v13_intermediates`）。

### 5.5 可复现性：确定性与重放

| 论断 | 测量 |
|---|---|
| 同一输入跑两次 → 科学内容相同 | 44 个可比文件，**0 个不同**（时间戳、运行时长、每次运行的路径已归一）（`v15_determinism`） |
| 标记排序与运行次序无关 | 两次运行排序一致（`v15_ranking`） |
| 记录的 `run_config.json` 可重放 | 不一致参数：**[]**（`v02_replay`） |
| 重放运行能复现产物 | 与记录运行的标记集一致（`v15_replay`） |
| 缓存能说明自己是根据什么构建的 | 19 个缓存基因树，19 个输入指纹；重跑是复用而不是重建（`v14_tree_cache`） |
| 随机种子被钉住而非交给时钟 | 两次运行的基因树与 support 值一致（`v15_seeds`、`v09_ufboot`） |

其中两条是新增的。重放一份快照过去是「名为重放的空操作」：`run_config.json` 按模块
嵌套参数并改名了其中几个字段，因此加载器会把 marker 目录、taxonomy 表、风险带当作
「未识别键」丢掉；路径按输出目录相对化写入、却被按当前工作目录解析，于是重放会去够
`../../../tmp/…`。而逐标记基因树缓存只按 marker id 复用，于是**别的序列**推出来的树可以
回答下一次运行的请求。现在：快照记录解析后的绝对路径，读取时展平回命令行词汇，并且
恢复「用户请求的」预算而不是自适应后的值；缓存条目携带其输入序列的 SHA-256。

### 5.6 稳健性：拒绝必须具体且大声

退化输入是被检验的，不是被绕开的：2-tip 与 1-tip 的 marker 文件、3 基因组的集合（低于
建树所需的 4-tip 下限）、全是不可用 profile 的 marker 目录、空的 `.faa`、非 FASTA 文件、
首行会被误当成数据的 taxonomy 表、非法参照树、缺失的 HMM 目录。每个用例都要求：非零退出
码、消息点名越界的文件、不出现 traceback、不产出空的占位产物
（`v16_degenerate.exit_code = 3`、`v16_three_tips.exit_code = 3`）。

通用的拒绝文案过去在 Phase 1 什么都没选中时仍然声称「所有标记都被判为 Level 3」；现在
两种情形分别说明，因为针对一种的建议对另一种是错的。

### 5.7 对植入阳性结果的检出（12 基因组基准）

成对设计：12 基因组 × 20 marker 的干净目录，以及同一目录植入 3 次跨目转移的版本。两次
运行都由套件真实执行；下列数字读自产物，而非先验断言
（`validation/results/chimera_metrics.json`）。

| 量 | 值 |
|---|---|
| 两次运行都被评分的 marker | 20 |
| 嵌合运行中被剔除（risk = 1.0）的植入 marker | **3 / 3** |
| 相对干净运行风险**上升**的植入 marker | 2 |
| 干净运行里就已达到最大风险的植入 marker | 1（`TIGR00061`）——这是被实测到的假阳性，按原样记录 |
| 成对运行之间发生变化的非植入 marker | 0 |
| 干净数据上的背景剔除率 | 0.1 |
| 剔除决策处的 precision / recall / F1 | 1.0 / 0.667 / 0.8 |

这个区分很重要：干净运行已经标记为阳性的 marker 不可能再「上升」，因此要求三个全部上升
要么不可满足，要么会把基线错误率藏起来。所以该用例断言的是：每个植入 marker 都被剔除、
旁观者保持不动、基线没有饱和——并把基线数字公布出来。`xfail` 只用于套件明确知道自己还
不能主张的结论；本次运行没有任何一处靠预期失败装成绿色。

## 6. 本次验收发现并修复的缺陷

| # | 症状（用户会看到什么） | 根因 | 钉住它的测试 |
|---|---|---|---|
| 1 | `--max-markers` 被接受、被记录、被忽略 | 自适应预算覆盖了用户值 | `v07_budget`、单元接线测试 |
| 2 | must-pass 打印「所有必需关系成立」，而它从没读过那棵树 | 把基因树**路径**当 Newick 解析 → 测试组数为 0 | `v11_gate_pass/fail` |
| 3 | 断言失败退出 3（运行时错误），与文档的 4 矛盾 | 聚合路径未映射该异常类型 | `v11_bypass`、退出码用例 |
| 4 | 所有基于 taxonomy 表的 marker 都评成 UNKNOWN | 跨阶元标志把 scope 推导的阶元也算作偏离 | `v04_format_b`、单元阶元权威性测试 |
| 5 | `--monophyly-rank` 被自动阶元悄悄替换 | scope 规则压过了文档承诺 | 单元 + `v04_ranks` |
| 6 | 7 对 1 的数据上一致 marker 全被判 risk 1.0 | 单系检验依赖定根；无信息量的 split 被计入 | `test_monophyly_split_semantics.py`、`v06` 风险带用例 |
| 7 | `--taxonomy-levels` 不起作用，kingdom 阶元被丢弃并告警 | 值被记录但从未传给表解析器 | `v04_custom_levels`、单元 parser 测试 |
| 8 | `--mol-type DNA` 接受了蛋白序列 | 该开关从未到达字母表校验 | `test_inert_switches_wired.py`、`v05_moltype` |
| 9 | 带 `--force` 的重跑可能复用由别的序列构建的基因树 | 缓存仅按 marker id 键控 | `v14_tree_cache`、单元 must-fail 对照 |
| 10 | `--config <上次运行>` 实际用了默认值 | 加载器看不懂嵌套的快照块 | `v02_replay`、`test_config_loader.py` |
| 11 | 重放会去够 `../../../tmp/…` 并在只读 `/tmp` 上崩溃 | 快照按输出目录相对化、读取时按 cwd 解析 | `v02_replay`、`v15_replay` |
| 12 | 每次运行看起来都不确定（`run_config.json`） | 归一化器的时间戳正则先破坏了 JSON 再解析 | `v15_determinism` |
| 13 | `-o` 指向已存在文件时抛裸 traceback | 未检查 `Path.is_file` | 退出码 2 用例 |
| 14 | 空 `.faa` 被当成一个基因组；不可读文件被静默忽略 | 加载时不报告种类与大小 | `v16_robustness_edges` |
| 15 | taxonomy 表的表头行变成一个幽灵 tip，进入每一次 tip 集比较 | 无表头嗅探 | `v04_*`、单元测试 |
| 16 | banner 把 pandas/numpy/scipy/rich 声称为本项目依赖 | 硬编码清单 | 改由发行版元数据派生 |
| 17 | 废弃别名被静默接受；`--tree` 的弃用提示永远不可能触发 | 与正主共用同一个 dest | `v05`/`v06`/`v07` 别名用例 |
| 18 | 退出 3 的文案断言了一个并未发生的原因 | 「什么都没选中」与「全被剔除」共用一条消息 | `v16_*`、`v07` 预算用例 |
| 19 | `PipelineConfig(output_dir=X)` 仍把报告写进 `./output` | 各 section 数据类各自带着 `output_dir`/`tmp_dir` 默认值，一次运行对“写到哪”有两个答案 | `test_config.py::TestOneAnswerPerRunLocation` |
| 20 | 暂存目录缺失被报成“hmmsearch 不在 PATH 上” | 输出路径的 `FileNotFoundError` 与工具缺失被混为一谈，结论指向了错的原因 | `test_run_scan_hmm_mode`；HMM 路径现在先建好自己的暂存目录 |
| 21 | 测试数据：taxonomy 错误、Newick 非法、嵌合目录与干净目录不配对、类别映射全 `NA` | 手工填值、写入器缺陷、缓存过期 | `01`–`03` 的构造期校验、`v17` 真值用例 |

## 7. 本次验收**没有**建立什么

* **与 GTDB-Tk 的等价性。** marker 输入是对 GTDB-Tk *布局*的重构（一 marker 一 FASTA、
  未比对序列、文件名即 marker id），由 hmmsearch 扫描这些本地拼装的 profile 得到；未复现 GTDB-Tk
  自身结合分类学的最佳命中选择。上面的数字度量的是 MarkerFinder 在布局一致输入上的
  判定，不是 GTDB-Tk 的标记选择。
* **流行病学意义上的检出率。** 12 基因组上 3 次跨目植入，外加实测背景率，是答案已知的
  阳性对照，不是真实 HGT 的样本。更大的估计来自 `tests/benchmark/`（30 基因组、GTDB
  r53 真值），因下载体积而由用户触发。
* **COG 生物学。** 功能类别列取自这些 profile 自身的注释行，因此它检验的是
  `--cog-category-map` 的通路。
* **占居率的极端分布。** 这套 profile 是保守标记集，在这组亲缘较近的基因组上最稀疏的天然
  marker 占居率仍为 0.875；更稀疏的输入只能靠截取模拟，并在
  `validation/data/markers/degenerate/README` 中如实标注。
* **Python 3.13 及以上。** ete3 在包导入时会 `import cgi`，而 CPython 3.13 已移除 `cgi`
  （PEP 594），此前导致 ete3 —— 以及所有 MAD / 单系性测量 —— 在更新的解释器上不可用。
  `markerfinder._cgi_compat` 现在只在标准库 `cgi` 确实缺失时安装最小替身，因此 3.13 与
  3.14 已纳入受支持范围，CI 矩阵逐个实跑。「不设上界」这一点本身在被测，而不是靠文档
  声明。残留风险：替身只覆盖 ete3 真正引用的那一个符号（`cgi.FieldStorage`），一旦被真正
  调用就抛 `NotImplementedError`；MarkerFinder 不使用 ete3 的 web 插件，该路径不会被触发。
* **可写的 `/tmp`。** MAFFT 用 `mktemp -dt` 建暂存目录，失败后并不会停止，而是带着空路径
  继续、产不出比对——流水线于是把结果报成「没有 marker 可比对」。套件把
  `TMPDIR`/`MAFFT_TMPDIR` 钉在 `validation/.work` 内，因此只读 `/tmp` 的容器不可能把一个
  环境问题变成错误的科学结论。

## 8. 如何复现本报告

```bash
conda env create -f environment.yml && conda activate markerfinder
# 先拼装 HMM profile（第三方模型，不入库）：python scripts/fetch_marker_db.py --auto
# 可选：从 NCBI 重建测试数据（需要网络）
python validation/run_validation.py --prepare
# 验收层，含慢速检出基准
python validation/run_validation.py --all -n 8
# 快速层（本报告同样引用其中的数字）
pytest tests
```

`run_validation.py` 会写出 `validation/results/{junit.xml,metrics.json,
capability_matrix.tsv,environment.txt}`；本报告引用这些文件，一致性测试会重新读取它们。
