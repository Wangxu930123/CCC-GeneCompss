# CCC-GeneCompass

## 简介

这是基于单细胞基础模型 [GeneCompass](https://github.com/xCompass-AI/GeneCompass) 的细胞互作分析工具。它利用深度学习和大模型方法来分析和研究细胞间的相互作用和通信。

本项目整合了：
1. 基于 CellChat 的细胞互作分析
2. 基于 CellPhoneDB 的细胞互作分析
3. 构建细胞互作分析金标准
4. 基于 GeneCompass 生成 Embeddings
5. 基于 GeneCompass 的细胞互作分析

## 主要功能

### 1. 细胞互作分析工具集成

本项目整合了两种主流的细胞互作分析工具：
- **CellChat**: 基于配体-受体相互作用的细胞通讯分析工具
- **CellPhoneDB**: 基于统计显著性检验的细胞互作分析工具

### 2. 金标准构建

通过整合 CellChat 和 CellPhoneDB 的分析结果，构建细胞互作金标准：
- 使用共识分数 (Consensus Score) 评估互作置信度
- 支持多种阈值方法 (分位数/绝对阈值)
- 生成完整的机器学习数据集

### 3. GeneCompass 模型应用

- 基于预训练的 GeneCompass 模型生成细胞/基因嵌入
- **使用基于 GeneCompass 的回归任务微调进行细胞互作预测**
  - 直接在 GeneCompass 模型上进行微调
  - 采用回归任务预测细胞互作强度分数
  - 保留 GeneCompass 模型的深层特征学习能力

## 项目结构

```
CCC-GeneCompass/
├── CellChatAnalysis/              # CellChat 分析工具
│   ├── h5ad_to_csv_fixed.py       # h5ad 转换为 CSV
│   ├── csv_to_rds_fixed.R         # CSV 转换为 RDS
│   ├── CellChatAnalysis_fixed.R   # CellChat 分析主脚本
│   ├── README.md                  # CellChat 使用文档
│   └── ...
├── CellPhoneAnalysis/             # CellPhoneDB 分析工具
│   ├── prepare_DEGs_h5ad.py       # 生成差异表达基因
│   ├── prepare_microenvs_h5ad.py  # 生成微环境文件
│   ├── CellPhoneAnalysis.py       # CellPhoneDB 分析主脚本
│   └── ...
├── genecompass/                   # GeneCompass 模型相关
│   ├── modeling_bert.py           # BERT 模型定义
│   ├── data_collator.py          # 数据整理器
│   └── ...
├── preprocess/                    # 数据预处理
│   ├── filter.py                  # 数据过滤
│   └── normalized.py             # 数据归一化
├── building_gold_standard_database.py  # 金标准构建
├── cell_cell_interaction.py       # 细胞互作分析
├── generate_embeddings.py         # 生成嵌入
└── requirements.txt               # 依赖包
```

## 快速开始

### 环境要求

- Python: 3.8+
- R: 4.0+
- 内存: 至少 16GB (推荐 32GB+)
- GPU: 可选，用于加速 GeneCompass 模型训练

### 安装依赖

```bash
# Python 环境
pip install -r requirements.txt

# R 环境 (用于 CellChat 和 CellPhoneDB)
# 参见各工具的 README 文档
```

## 完整分析流程

### 步骤 1: 数据预处理

```bash
# 数据过滤
python preprocess/filter.py --input data.h5ad --output filtered.h5ad

# 数据归一化
python preprocess/normalized.py --input filtered.h5ad --output normalized.h5ad
```

### 步骤 2: 运行 CellChat 分析

```bash
cd CellChatAnalysis

# 转换 h5ad 为 CSV
python h5ad_to_csv_fixed.py \
  --input ../normalized.h5ad \
  --output ./cellchat_output \
  --celltype_col cell_type

# 创建 Seurat 对象
Rscript csv_to_rds_fixed.R ./cellchat_output

# 运行 CellChat 分析
Rscript CellChatAnalysis_fixed.R \
  ./cellchat_output/seurat_obj.rds \
  ./cellchat_output \
  cell_type
```

### 步骤 3: 运行 CellPhoneDB 分析

```bash
cd CellPhoneAnalysis

# 生成差异表达基因文件
python prepare_DEGs_h5ad.py \
  --h5ad ../normalized.h5ad \
  --outdir ./cpdb_output \
  --groupby cell_type

# 生成微环境文件
python prepare_microenvs_h5ad.py \
  --h5ad ../normalized.h5ad \
  --output ./cpdb_output/microenv.tsv \
  --groupby cell_type

# 运行 CellPhoneDB 分析
python CellPhoneAnalysis.py \
  --h5ad ../normalized.h5ad \
  --cpdb ./v5.0.0/cellphonedb.zip \
  --degs ./cpdb_output/DEGs.tsv \
  --microenv ./cpdb_output/microenv.tsv \
  --groupby cell_type \
  --outdir ./cpdb_output
```

### 步骤 4: 构建金标准

```bash
cd ..

python building_gold_standard_database.py
```

此脚本将整合 CellChat 和 CellPhoneDB 的结果，生成金标准数据集。

### 步骤 5: 生成 Embeddings

```bash
python generate_embeddings.py \
  --dataset_path ./normalized_data/ \
  --model_path ./pretrained_models/GeneCompass_Base \
  --token_dict_path ./prior_knowledge/human_mouse_tokens.pickle \
  --output_path ./embeddings/embeddings.pickle
```

### 步骤 6: 细胞互作分析

```bash
python cell_cell_interaction.py \
  --dataset_path ./normalized_data/ \
  --gold_standard_path ./gold_standard/machine_learning_dataset.csv \
  --token_dict_path ./prior_knowledge/human_mouse_tokens.pickle \
  --model_path ./pretrained_models/GeneCompass_Base \
  --embeddings_path ./embeddings/embeddings.pickle \
  --output_dir ./outputs/
```

## 输出文件

### CellChat 输出
- `cellchat_communication.csv` - 细胞通讯结果
- `cellchat_pathways.csv` - 信号通路结果
- `cell_interaction_strength_matrix.csv` - 互作强度矩阵
- `communication_network.png` - 通讯网络图

### CellPhoneDB 输出
- `significant_means.txt` - 显著互作结果
- `means.txt` - 所有互作均值
- `pvalues.txt` - 显著性检验结果

### 金标准输出
- `complete_labeled_interactions.csv` - 完整标记矩阵
- `machine_learning_dataset.csv` - 机器学习数据集
- `gold_standard_interactions.csv` - 金标准正样本
- `dataset_statistics.csv` - 数据集统计信息

### GeneCompass 输出
- `embeddings.pickle` - 细胞/基因嵌入
- `trained_model/` - 训练好的模型
- `predictions.csv` - 预测结果
- `metrics.json` - 评估指标

### GeneCompass预训练模型获取

[GeneCompass](https://github.com/xCompass-AI/GeneCompass)预训练模型可以通过下面的链接获取：

将pretrained_model目录置于主路径下（`./pretrained_models/GeneCompass_Small`，`./pretrained_models/GeneCompass_Base`）

| Model             | Description                         | Download                                           |
| ----------------- | ----------------------------------- | -------------------------------------------------- |
| GeneCompass_Small | Pretrained on 6-layer GeneCompass.  | [Link](https://www.scidb.cn/en/anonymous/SUZOdk1y) |
| GeneCompass_Base  | Pretrained on 12-layer GeneCompass. | [Link](https://www.scidb.cn/en/anonymous/SUZOdk1y) |

## 引用

如果您在研究中使用了本工具，请引用以下论文：

- GeneCompass:Yang, X., Liu, G., Feng, G. *et al.* GeneCompass: deciphering universal gene regulatory mechanisms with a knowledge-informed cross-species foundation model. *Cell Res* **34**, 830–845 (2024). https://doi.org/10.1038/s41422-024-01034-y
- CellChat: Jin, S., et al. (2021). Inferring cell-cell communication by integrating ligand-receptor, signaling gene, and TF-target networks. Nature Protocols.
- CellPhoneDB: Vento-Tormo, R., et al. (2018). Single-cell reconstruction of developmental trajectories during human endometriosis. Science.

## 许可证

本项目遵循原项目 GeneCompass 的许可证。

## 联系方式

如有问题或建议，请提交 Issue 或 Pull Request。
