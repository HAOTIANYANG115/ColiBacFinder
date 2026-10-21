# ColiBacFinder 0.1.0

大肠杆菌细菌素相关基因组候选区域筛查工具，研究原型版。

输入一个核酸 FASTA 文件或一个目录，无需先提供 FAA/GFF。使用固定参考蛋白对核酸序列进行 TBLASTN 搜索，自动生成样本总表、候选位置、备选命中、系统组分和邻域序列。只读取输入目录第一层 `.fa/.fna/.fasta`（支持 `.gz`），忽略其他文件。

**结果是序列与邻域证据，不是细菌素活性鉴定。** “near_reference”不是经过独立数据校准的阳性标准；免疫单独命中、共享转运蛋白、部分结构域命中不能当作完整细菌素系统。标签数量不是独立基因数量。

## 安装

当前可通过 Conda 安装依赖，再从源码安装 ColiBacFinder；尚未发布独立 Conda 软件包。

```bash
git clone https://github.com/HAOTIANYANG115/ColiBacFinder.git
cd ColiBacFinder
conda env create -f environment.yml
conda activate colibacfinder
python -m pip install .
colibacfinder --version
colibacfinder db-check
colibacfinder -h
colibacfinder run -h
```

私有仓库需要已获授权的 GitHub 账号才能克隆。

依赖：Python ≥3.9，NCBI BLAST+（`makeblastdb`、`tblastn` 必须在 PATH）。本工具没有第三方 Python 运行依赖。附带冻结数据库，无需预先运行 Prokka、ABRicate 或提供 FAA/GFF。

## 一条命令分析

```bash
colibacfinder run -i /path/to/genome_fasta -o /path/to/new_results -t 8
```

`-i`：核酸 FASTA 文件或目录；`-o`：新结果目录；`-t`：每个样本的线程数。样本按顺序处理。

后台运行：

```bash
nohup colibacfinder run -i /path/to/genome_fasta -o /path/to/new_results -t 8 > colibacfinder.log 2>&1 &
tail -f colibacfinder.log
```

输出目录必须不存在，避免覆盖旧结果。退出码：0 全部完成；1 存在失败样本；2 参数/准备错误；130 中断。当前版本不支持断点续跑，重跑请用新目录。不会修改原始输入。

## 主要输出

| 文件 | 内容 |
|---|---|
| `01_sample_summary.tsv` | 每个样本一行：候选区域数、近参考本体/免疫标签、系统候选、失败信息 |
| `02_all_HSPs.tsv` | 全部原始搜索命中的分级、坐标、参考登录号、参考等级 |
| `03_candidate_loci.tsv` | 合并后的候选位置；最佳与次佳参考、得分差、歧义、边界标记 |
| `04_locus_evidence.tsv` | 各位置所有支持命中，保留备选标签 |
| `05_proximity_pairs.tsv` | 本体/免疫候选的距离、重叠、方向；仅表示邻近 |
| `06_system_candidates.tsv` | 按系统列出候选与组分支持状态 |
| `07_system_components.tsv` | 组分逐条证据；弱命中不能补足近参考组分 |
| `08_label_frequency.tsv` | 标签频率，以成功样本为分母，类别可重叠 |
| `samples/S00001/` 等 | 样本原始命中、日志、contig 映射、候选邻域序列及坐标 |
| `input_manifest.tsv` | 样本名、输入路径/校验和与工作目录对应表 |
| `parameters.json` | 软件、数据库、搜索工具版本、阈值和来源信息 |
| `reference_snapshot/` | 本次实际使用的数据库快照，含来源及文献字段 |

邻域默认两侧各 5 kb；间距默认 10 kb。Microcin V 的前体与免疫命中允许重叠，单独记录重叠长度。候选数量不等于预测完整蛋白数量；`coding_integrity=not_assessed`。

## 固定默认规则

搜索 E-value ≤1e-3；入筛 E-value ≤1e-5。按优先顺序：

| 分级 | identity | 参考覆盖 | 最短比对 |
|---|---:|---:|---:|
| near_reference | ≥80% | ≥80% | 不额外限制 |
| remote | ≥35% | ≥70% | 40 aa |
| partial | ≥35% | ≥10% | 60 aa |

覆盖为单个 HSP 的参考坐标跨度，不是多 HSP 合并覆盖。当前阈值沿用探索流程，尚未做独立性能校准。不同角色和参考等级始终保留。相同位置存在多个类型标签时保留歧义，不强行选第一条。

Microcin V、J25、B17 的 `all_model_components_near_reference` 仅代表所有模型组分都有近参考邻域证据，未验证 ORF 完整性、表达或活性。H47/M/I47 及 colicin 模型有不同覆盖范围，详见输出 `system_models.json`。转运/成熟组分单列，不能增加本体数量。

## 重出报告与独立评估

```bash
colibacfinder report --from-run /path/to/old_run -o /path/to/new_report
colibacfinder evaluate --run /path/to/run --truth truth.tsv -o /path/to/evaluation
```

重出报告不重跑 BLAST，核对输入校验和与数据库指纹；修改阈值可用 `--identity`、`--coverage` 等。旧版21号脚本输出可尝试读取，但没有原始输入校验和，报告会标记这一限制。新的规则更明确地处理辅助组分，因此不保证计数与旧脚本完全一致。

`truth.tsv` 必须含 `sample,target,truth,evidence` 四列（制表符），truth 为 0/1。支持 `ANY_CORE`、`ANY_NONEXPORT_IMMUNITY`、`system:Microcin_V`、`system:Microcin_J25`、`system:Microcin_B17`。证据需对应同一预测终点；不应把“培养未出现抑菌圈”直接当作“基因不存在”。未知状态不填为阴性。失败/缺失样本自动排除，并保留排除记录。

独立评估前先冻结数据库和阈值，留出未参与建库与调参的阳性、阴性样本，分层评估本体/免疫/系统，报告 TP/FP/TN/FN、敏感度、特异度和精确率。自匹配与模拟序列只用于程序自检，不能用于报告真实准确率。

## 环状序列

仅对已确认闭合序列可提供 `--circular-map circular.tsv`，两列 `sample`、`contig`。它用于跨原点邻域距离与序列提取，不会自动拼接跨原点的 TBLASTN 命中。不能仅根据质粒名称认定闭合。

## 版本与范围

固定库：`bacteriocin_unified_v0.1_20261020_194442_905507`，84 条独立参考蛋白、92 条注释关系。包括已审阅主参考、候选、片段和辅助参考，不能全部称作已验证参考。原始元数据和证据链接原样保留；可在 `reference_snapshot/metadata/reference_links.tsv` 查登录号和证据，在 `source_metadata/` 查详细来源。

软件不进行物种鉴定、质粒定位鉴定、ORF 完整性确认、毒素活性判断或精确亚型确认；不包含独立 nucleotide BLAST 路径。碎片化、缺失、短肽统计敏感性和不完整参考覆盖可能造成漏检。科研原型，发布许可及作者信息须由项目所有者在公开发布前确定；未擅自改变原始数据授权。
