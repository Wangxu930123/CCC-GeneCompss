# CCC-GeneCompass 快速入门指南

## 简介

CCC-GeneCompass 是一个基于 GeneCompass 单细胞基础模型的细胞互作分析工具。本指南将帮助您快速开始使用该工具进行细胞互作分析。

## 分析流程概述

```
单细胞数据 (.h5ad)
    │
    ├─────────────┬─────────────┬─────────────┐
    ▼             ▼             ▼             ▼
数据预处理    CellChat分析   CellPhoneDB分析
    │             │             │
    ▼             ▼             ▼
归一化数据    互作结果      互作结果
    │             │             │
    └─────────────┴─────────────┘
                  │
                  ▼
            构建金标准
                  │
                  ▼
          生成 Embeddings
                  │
                  ▼
        GeneCompass 训练与预测
                  │
                  ▼
            细胞互作网络
```

## 环境准备

### 1. Python 环境

```bash
# 创建 conda 环境
conda create -n ccc python=3.10
conda activate ccc

# 安装依赖
pip install -r requirements.txt
```

### 2. R 环境 (用于 CellChat 和 CellPhoneDB)

```bash
# 安装必要的 R 包
R

# 在 R 中执行
install.packages(c("Seurat", "dplyr", "ggplot2", "patchwork"))
if (!requireNamespace("BiocManager", quietly = TRUE))
    install.packages("BiocManager")
BiocManager::install("ComplexHeatmap")

# 安装 CellChat
install.packages("devtools")
devtools::install_github("sqjin/CellChat")

# 安装 CellPhoneDB
# 下载 CellPhoneDB 并安装到 conda 环境
# 详见 CellPhoneDB 官方文档
```

## 数据准备

### 1. 下载预训练模型

从 GeneCompass 官方仓库下载预训练模型：

- GeneCompass_Base: [下载链接](https://www.scidb.cn/en/anonymous/SUZOdk1y)
- GeneCompass_Small: [下载链接](https://www.scidb.cn/en/anonymous/SUZOdk1y)

将模型文件放置在 `./pretrained_models/` 目录下。

### 2. 准备先验知识

从 GeneCompass 仓库下载以下文件，放置在 `./prior_knowledge/` 目录下：

- `human_mouse_tokens.pickle`
- `human_gene_median_after_filter.pickle`
- 其他必要的知识文件

### 3. 准备单细胞数据

确保您的单细胞数据是 `.h5ad` 格式，并包含：
- 表达矩阵 (`adata.X`)
- 细胞元数据 (`adata.obs`)，包含细胞类型列
- 基因信息 (`adata.var`)

## 完整分析流程

### 步骤 1: 数据预处理 (可选)

如果您的数据已经归一化，可以跳过此步骤。

```bash
# 数据过滤
python preprocess/filter.py

# 数据归一化
python preprocess/normalized.py
```

### 步骤 2: 运行 CellChat 分析

```bash
cd CellChatAnalysis

# 2.1 转换 h5ad 为 CSV
python h5ad_to_csv_fixed.py \
  --input ../your_data.h5ad \
  --output ./cellchat_output \
  --celltype_col cell_type

# 2.2 创建 Seurat 对象
Rscript csv_to_rds_fixed.R ./cellchat_output 3 200

# 2.3 运行 CellChat 分析
Rscript CellChatAnalysis_fixed.R \
  ./cellchat_output/seurat_obj.rds \
  ./cellchat_output \
  cell_type \
  20 \
  4 \
  10 \
  FALSE
```

**输出文件**：
- `cellchat_communication.csv` - 细胞通讯结果
- `cellchat_pathways.csv` - 信号通路结果
- `cell_interaction_strength_matrix.csv` - 互作强度矩阵

### 步骤 3: 运行 CellPhoneDB 分析

```bash
cd CellPhoneAnalysis

# 3.1 生成差异表达基因文件
python prepare_DEGs_h5ad.py \
  --mode degs_only \
  --h5ad ../your_data.h5ad \
  --outdir ./cpdb_output \
  --groupby cell_type \
  --n_top_genes 250

# 3.2 生成微环境文件
python prepare_microenvs_h5ad.py \
  --h5ad ../your_data.h5ad \
  --output ./cpdb_output/microenv.tsv \
  --groupby cell_type

# 3.3 运行 CellPhoneDB 分析
python CellPhoneAnalysis.py \
  --h5ad ../your_data.h5ad \
  --cpdb ./v5.0.0/cellphonedb.zip \
  --degs ./cpdb_output/DEGs.tsv \
  --microenv ./cpdb_output/microenv.tsv \
  --groupby cell_type \
  --outdir ./cpdb_output
```

**输出文件**：
- `significant_means.txt` - 显著互作结果
- `means.txt` - 所有互作均值
- `pvalues.txt` - 显著性检验结果

### 步骤 4: 构建金标准

```bash
cd ..

python building_gold_standard_database.py \
  --cellchat_dir ./CellChatAnalysis/cellchat_output/ \
  --cpdb_dir ./CellPhoneAnalysis/cpdb_output/ \
  --output_dir ./gold_standard/ \
  --threshold_method quantile \
  --threshold_value 0.7
```

**输出文件**：
- `complete_labeled_interactions.csv` - 完整标记矩阵
- `machine_learning_dataset.csv` - 机器学习数据集
- `gold_standard_interactions.csv` - 金标准正样本
- `dataset_statistics.csv` - 数据集统计信息
- 多种可视化图表 (`.png`)

### 步骤 5: 生成 Embeddings

```bash
python generate_embeddings.py \
  --dataset_path ./normalized_data/TabulaSapiens/tabula_sapiens_liver/ \
  --model_path ./pretrained_models/GeneCompass_Base \
  --token_dict_path ./prior_knowledge/human_mouse_tokens.pickle \
  --output_path ./embeddings/embeddings.pickle \
  --batch_size 128 \
  --gpu_ids "0"
```

**注意**：此步骤可能需要较长时间，取决于数据集大小。

### 步骤 6: 细胞互作分析

```bash
python cell_cell_interaction.py
```

**输出文件**：
- `best_model.pt` - 训练好的最佳模型
- `evaluation_results.pkl` - 评估结果
- `predictions.csv` - 预测结果（如果添加推理代码）

## 输出文件说明

### CellChat 输出
- `cellchat_communication.csv`: 每行一个细胞通讯事件，包含发送者、接收者、互作名称和概率
- `cellchat_pathways.csv`: 信号通路及其通讯概率
- `cell_interaction_strength_matrix.csv`: 细胞类型之间的互作强度矩阵
- `communication_network.png`: 细胞通讯网络可视化
- `communication_heatmap.png`: 互作强度热图

### CellPhoneDB 输出
- `significant_means.txt`: 显著的配体-受体互作及强度
- `means.txt`: 所有配体-受体互作的平均强度
- `pvalues.txt`: 显著性检验的 p 值

### 金标准输出
- `complete_labeled_interactions.csv`: 包含所有细胞对及其标签的完整数据集
- `machine_learning_dataset.csv`: 用于机器学习的数据集，包含特征和标签
- `gold_standard_interactions.csv`: 仅包含高置信度的金标准正样本
- `dataset_statistics.csv`: 数据集的统计信息
- `complete_interaction_matrix.png`: 完整互作矩阵热图
- `gold_standard_labels.png`: 金标准标签热图
- `score_distributions.png`: 共识分数分布图
- `gold_standard_network.png`: 金标准互作网络图

### GeneCompass 输出
- `embeddings.pickle`: 细胞/基因的嵌入向量
- `best_model.pt`: 训练好的最佳模型
- `evaluation_results.pkl`: 评估指标和结果

## 常见问题

### Q1: 内存不足怎么办？

**A**:
- 减小 `generate_embeddings.py` 中的 `--batch_size` 参数
- 减小 CellChat 的并行 worker 数量
- 在 CellPhoneDB 分析中使用数据降采样

### Q2: 基因名格式不匹配？

**A**:
- 检查数据使用的是 Ensembl ID 还是 Gene Symbol
- CellChat 和 CellPhoneDB 通常使用 Gene Symbol
- 如果是 Ensembl ID，需要先转换为 Gene Symbol

### Q3: 没有检测到任何细胞互作？

**A**:
- 检查细胞类型数量（至少需要 2 种）
- 检查每种细胞类型的细胞数是否足够
- 降低 `min_cells` 参数
- 检查基因名格式是否与数据库匹配

### Q4: 训练时间过长？

**A**:
- 使用 GPU 加速
- 减小数据集规模
- 调整模型复杂度
- 使用 GeneCompass_Small 模型而不是 GeneCompass_Base

### Q5: 如何评估模型性能？

**A**:
- 查看训练和验证集的准确率
- 使用评估指标：准确率、精确率、召回率、F1 分数
- 查看混淆矩阵
- 与传统方法（CellChat、CellPhoneDB）的结果比较

## 下一步

- 查看详细的 `README.md` 了解更多功能
- 查看各子目录的 `README.md` 了解详细使用方法
- 根据您的具体需求调整参数

## 引用

如果您在研究中使用了本工具，请引用：

- GeneCompass: [论文链接]
- CellChat: Jin, S., et al. (2021). Nature Protocols.
- CellPhoneDB: Vento-Tormo, R., et al. (2018). Science.

## 联系方式

如有问题或建议，请提交 Issue 或 Pull Request。
