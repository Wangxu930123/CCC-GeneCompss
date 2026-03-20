# CellChat细胞互作分析工具（修复版）

## 概述

本目录包含CellChat细胞互作分析所需的完整流程工具，从h5ad文件转换到最终的CellChat分析结果。

## 文件说明

### 主要脚本

| 文件 | 功能 | 说明 |
|------|------|------|
| `verify_gene_format.py` | 验证基因格式 | 检查h5ad中的基因格式，判断是否需要转换 |
| `generate_gene_mapping.py` | 生成基因映射 | 从h5ad提取Ensembl ID并查询Gene Symbol，生成映射文件 |
| `convert_ensembl_to_symbol.py` | 转换基因名 | 使用映射文件将h5ad中的Ensembl ID转换为Gene Symbol |
| `run_cellchat_analysis.bat/sh` | 自动化脚本 | 一键完成整个CellChat分析流程（包括基因转换） |
| `h5ad_to_csv_fixed.py` | h5ad转稀疏矩阵 | 将h5ad文件转换为CellChat分析所需的格式 |
| `csv_to_rds_fixed.R` | 稀疏矩阵转Seurat对象 | 将稀疏矩阵和元数据创建为Seurat RDS文件 |
| `CellChatAnalysis_fixed.R` | CellChat分析 | 从Seurat对象进行CellChat细胞互作分析 |

### 原始脚本（已修复）

- `h5ad_to_csv.py` - 原始版本（已废弃，请使用`_fixed`版本）
- `csv_to_rds.R` - 原始版本（已废弃，请使用`_fixed`版本）
- `CellChatAnalysis.R` - 原始版本（已废弃，请使用`_fixed`版本）

## 修复的问题

### 1. 硬编码路径 ✅
- **问题**: 所有路径都是Linux服务器路径
- **修复**: 添加命令行参数支持，可在任何系统运行

### 2. 基因名处理不一致 ✅
- **问题**: 基因名处理逻辑混乱，可能丢失信息
- **修复**: 统一基因标识符处理，支持Ensembl和Gene Symbol自动检测

### 3. 线粒体基因检测失效 ✅
- **问题**: 基因名大小写转换后模式匹配失效
- **修复**: 根据基因标识符类型选择合适的线粒体基因模式

### 4. 并行计算不兼容Windows ✅
- **问题**: 硬编码`"multicore"`，Windows不支持
- **修复**: 自动检测操作系统，Windows使用`"multisession"`

### 5. 数据库过滤过于严格 ✅
- **问题**: 排除所有含`-`、`_`的基因，丢失大量有效配体-受体对
- **修复**: 只排除明显的复合物（含COMPLEX字样），可通过参数控制

### 6. 错误处理过于宽松 ✅
- **问题**: 过多tryCatch掩盖真正问题
- **修复**: 记录详细错误信息，改进参数fallback逻辑

### 7. 输出格式不兼容 ✅
- **问题**: CSV格式与CellPhoneDB不一致
- **修复**: 添加CellPhoneDB兼容的配体-受体对输出

### 8. 缺少数据验证 ✅
- **问题**: 没有验证基因标识符格式、细胞类型一致性
- **修复**: 添加完整的数据验证和诊断信息

## 使用方法

### 重要提示：Ensembl基因ID处理

CellChat数据库使用**Gene Symbol**格式（如`TGFB1`, `TGFBR1`），如果您的数据使用**Ensembl ID**（如`ENSG00000000003`），需要先转换为Gene Symbol。

**检查数据格式**:
```bash
# 方法1: 使用验证工具（推荐）
python verify_gene_format.py your_data.h5ad

# 方法2: 手动检查
python -c "import scanpy as sc; adata=sc.read_h5ad('your_data.h5ad'); print(adata.var_names[:10])"
```

如果看到`ENSG00000000003`这样的ID，请使用以下两种方法之一：

#### 方法1: 自动转换（推荐）

如果Seurat对象中包含Gene Symbol信息，CellChat分析脚本会自动检测并转换。

#### 方法2: 使用映射文件（适用于纯Ensembl ID数据）

**步骤0: 生成Ensembl ID到Gene Symbol的映射文件**

```bash
# 安装mygene包
pip install mygene

# 生成映射文件（需要网络连接）
python generate_gene_mapping.py \
  --input your_data.h5ad \
  --output gene_mapping.tsv
```

**步骤0.5: 转换h5ad文件中的基因名**

```bash
python convert_ensembl_to_symbol.py \
  --input your_data.h5ad \
  --output your_data_with_symbols.h5ad \
  --mapping gene_mapping.tsv
```

然后在后续步骤中使用`your_data_with_symbols.h5ad`。

### 步骤1: h5ad转稀疏矩阵

```bash
python h5ad_to_csv_fixed.py \
  --input <path/to/input.h5ad> \
  --output <path/to/output_dir> \
  --celltype_col cell_type
```

**参数说明**:
- `--input`: 输入h5ad文件路径（必需）
- `--output`: 输出目录路径（必需）
- `--celltype_col`: 细胞类型列名（默认: cell_type）

**输出文件**:
- `sparse_matrix.mtx` - 稀疏矩阵
- `barcodes.tsv` - 细胞条形码
- `genes.tsv` - 基因名列表
- `metadata.csv` - 元数据
- `conversion_summary.txt` - 转换摘要

### 步骤2: 创建Seurat对象

```bash
Rscript csv_to_rds_fixed.R <output_dir> [min_cells] [min_features]
```

**参数说明**:
- `output_dir`: 输入文件所在目录（必需）
- `min_cells`: 每个基因的最小细胞数（默认: 3）
- `min_features`: 每个细胞的最小基因数（默认: 200）

**输出文件**:
- `seurat_obj.rds` - Seurat对象
- `QC_plot.png` - 质量控制图
- `seurat_conversion_summary.txt` - 转换摘要

### 步骤3: CellChat分析

```bash
Rscript CellChatAnalysis_fixed.R \
  <input_rds> \
  <output_dir> \
  [cell_type_col] \
  [min_cells] \
  [workers] \
  [nboot] \
  [keep_complex]
```

**参数说明**:
- `input_rds`: 输入Seurat RDS文件路径（必需）
- `output_dir`: 输出目录路径（必需）
- `cell_type_col`: 细胞类型列名（默认: cell_type）
- `min_cells`: 最小细胞数阈值（默认: 20）
- `workers`: 并行worker数量（默认: 4）
- `nboot`: bootstrap次数（默认: 10）
- `keep_complex`: 是否保留复合物配体-受体对（默认: FALSE）

**输出文件**:
- `cellchat_result.rds` - CellChat对象
- `cellchat_communication.csv` - 细胞通讯结果
- `cellchat_pathways.csv` - 信号通路结果
- `cell_type_info.csv` - 细胞类型信息
- `cell_interaction_strength_matrix.csv` - 细胞互作强度矩阵
- `cellchat_lr_pairs.csv` - **配体-受体对（CellPhoneDB兼容）**
- `communication_network.png` - 通讯网络图
- `communication_heatmap.png` - 通讯热力图
- `pathway_bubble.png` - 信号通路气泡图
- `ligand_receptor_network.png` - 配体-受体互作图
- `session_info.txt` - 会话信息

## 完整示例

### 情况1: 数据使用Gene Symbol（推荐）

```bash
# 1. 转换h5ad
python h5ad_to_csv_fixed.py \
  --input G:/DATA/SingleCell/tabula_sapiens_liver.h5ad \
  --output ./cellchat_output

# 2. 创建Seurat对象
Rscript csv_to_rds_fixed.R ./cellchat_output 3 200

# 3. 运行CellChat分析
Rscript CellChatAnalysis_fixed.R \
  ./cellchat_output/seurat_obj.rds \
  ./cellchat_output \
  cell_type \
  20 \
  4 \
  10 \
  FALSE
```

### 情况2: 数据使用Ensembl ID

```bash
# 0. 生成基因映射文件（需要网络连接）
pip install mygene
python generate_gene_mapping.py \
  --input your_data.h5ad \
  --output gene_mapping.tsv

# 0.5. 转换基因名
python convert_ensembl_to_symbol.py \
  --input your_data.h5ad \
  --output your_data_with_symbols.h5ad \
  --mapping gene_mapping.tsv

# 1. 转换h5ad
python h5ad_to_csv_fixed.py \
  --input your_data_with_symbols.h5ad \
  --output ./cellchat_output

# 2. 创建Seurat对象
Rscript csv_to_rds_fixed.R ./cellchat_output 3 200

# 3. 运行CellChat分析
Rscript CellChatAnalysis_fixed.R \
  ./cellchat_output/seurat_obj.rds \
  ./cellchat_output \
  cell_type \
  20 \
  4 \
  10 \
  FALSE
```

## 与CellPhoneDB集成

### CellChat配体-受体对输出格式

`cellchat_lr_pairs.csv`包含以下列：
- `ligand` - 配体基因
- `receptor` - 受体基因
- `source_db` - 数据源标识（"CellChat"）

该格式与CellPhoneDB分析结果兼容，可用于：
1. 构建CellPhoneDB + CellChat共识金标准
2. 比较两个数据库的预测结果
3. 合并分析结果

### 使用示例

```python
import pandas as pd

# 读取CellChat配体-受体对
cellchat_lr = pd.read_csv('./cellchat_output/CellChat_Results/cellchat_lr_pairs.csv')

# 读取CellPhoneDB配体-受体对
cpdb_lr = pd.read_csv('./cellphone_output/degs_analysis_relevant_interactions.tsv', sep='\t')

# 创建共识
cellchat_set = set(zip(cellchat_lr['ligand'], cellchat_lr['receptor']))
cpdb_set = set(zip(cpdb_lr['gene_name'], cpdb_lr['gene_name_1']))

# 交集（共识）
consensus = cellchat_set & cpdb_set

# 输出共识配体-受体对
consensus_df = pd.DataFrame(list(consensus), columns=['ligand', 'receptor'])
consensus_df['source'] = 'consensus'
consensus_df.to_csv('./consensus_lr_pairs.csv', index=False)
```

## 环境要求

### Python环境

```bash
pip install anndata scipy pandas numpy
```

### R环境

```r
install.packages(c("Seurat", "CellChat", "ggplot2", "patchwork",
                   "dplyr", "future", "ComplexHeatmap", "RColorBrewer", "stringr"))
```

或使用conda:

```bash
conda install -c bioconda r-seurat r-cellchat
```

## 常见问题

### Q1: Windows下R脚本报错"multicore not supported"

**A**: 使用`_fixed`版本的脚本，会自动检测操作系统并使用`"multisession"`。

### Q2: 基因名格式不匹配

**A**: 新脚本会自动检测Ensembl和Gene Symbol格式，并正确转换为大写。

### Q3: 内存不足

**A**:
- 减少`workers`参数（如设置为2或4）
- 减少`nboot`参数（如设置为5）
- 在csv_to_rds中增加`min_cells`和`min_features`过滤

### Q4: 细胞类型名称不一致

**A**: 使用`--celltype_col`或`cell_type_col`参数指定正确的列名。

### Q5: 没有检测到任何细胞互作

**A**:
- 检查细胞类型数量是否足够（至少2种）
- 检查每个细胞类型的细胞数是否达标（默认≥20）
- 降低`min_cells`参数
- 增加`keep_complex=TRUE`以保留更多配体-受体对
- **检查基因名格式是否与数据库匹配（CellChat数据库使用大写Gene Symbol）**
- 如果使用Ensembl ID，请参考"重要提示：Ensembl基因ID处理"部分进行转换

### Q7: 预处理后数据维度为0 x N（没有基因与数据库匹配）

**A**: 这是最常见的问题，原因通常是：
- 数据使用Ensembl ID（如`ENSG00000000003`），但CellChat数据库使用Gene Symbol（如`TGFB1`）
- 解决方法：使用`generate_gene_mapping.py`和`convert_ensembl_to_symbol.py`进行转换
- 参考README中的"情况2: 数据使用Ensembl ID"部分

### Q6: 预处理后数据维度不匹配错误

**A**: 这是CellChat的已知问题，修复版脚本已自动处理。如果仍有问题：
- 检查输入数据是否包含足够的细胞（建议>1000个）
- 检查细胞类型命名是否一致
- 尝试降低`min_cells`参数

## 作者

GeneCompass团队

## 版本历史

- v2.0 (2026-02-26) - 修复所有已知问题，添加命令行参数支持
- v1.0 - 初始版本
