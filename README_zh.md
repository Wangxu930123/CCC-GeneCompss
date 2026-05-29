<<<<<<< HEAD
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
=======
# CCC-GeneCompass：基于GeneCompass大模型的细胞互作分析

**利用大规模预训练语言模型进行单细胞转录组细胞间通讯预测**

[![Python 3.11+](https://img.shields.io/badge/python-3.11+-blue.svg)](https://www.python.org/)
[![PyTorch](https://img.shields.io/badge/PyTorch-2.5.0-red.svg)](https://pytorch.org/)
[![License](https://img.shields.io/badge/License-MIT-green.svg)](LICENSE)
[English](README.md)

---

## 概述

CCC-GeneCompass 基于 **GeneCompass**——一个在超过1亿个单细胞转录组上预训练的BERT架构大模型——预测细胞间相互作用（CCI）强度。分析流程包括：

1. 从单细胞数据计算**细胞类型聚合表达谱**
2. 基于 CellChat + CellPhoneDB v5 共识构建**金标准**
3. 微调 GeneCompass 对细胞类型相互作用强度进行排序
4. 通过**5折交叉验证**和**Spearman秩相关系数ρ**进行评估

## 核心特性

- **确定性**：细胞类型均值表达，消除单细胞随机采样噪声
- **科学严谨**：Spearman ρ + Bootstrap 95%置信区间 + 置换检验p值
- **多GPU支持**：DataParallel训练，最多支持4块GPU
- **端到端流程**：数据预处理 → 金标准构建 → 训练 → 交叉验证 → 可视化
- **模块化设计**：各步骤独立脚本，便于定制

## 安装

```bash
pip install -r requirements.txt

# CellChat 的 R 依赖
R -e 'install.packages(c("devtools","NMF","circlize","ComplexHeatmap"))'
R -e 'devtools::install_github("jinworks/CellChat")'

# CellPhoneDB v5 数据库
# 从 https://github.com/ventolab/CellphoneDB 下载 cellphonedb.zip
```

## 数据准备

### 1. 原始单细胞数据 (h5ad)
要求：`.obs` 中包含 `cell_type` 列，`.var` 中可通过 `feature_name` 或 `var_names` 获取基因符号。

### 2. 预训练模型
下载 GeneCompass_Base 模型文件，放入：
```
pretrained_models/
├── pytorch_model.bin   # ~1.1GB
└── config.json
```

### 3. 知识文件（用于token化）
```
prior_knowledge/
├── human_mouse_tokens.pickle
└── public/
    └── human_gene_median_after_filter.pickle
```

### 4. CellPhoneDB v5 数据库
```
CellPhoneAnalysis/v5.0.0/
├── cellphonedb.zip
├── gene_input.csv
└── protein_input.csv
```

## 分析流程

```bash
# ================================================================
# 单器官完整分析流程
# ================================================================
# 输入：  原始数据 h5ad（需包含 cell_type 列）
# 输出：  results/ 目录（5折交叉验证指标 + 模型 + 可视化）
# ================================================================

ORGAN=pancreas
RAW_H5AD=/path/to/original/${ORGAN}.h5ad       # 原始单细胞表达数据

# ====== Step 0: 数据预处理 ======
#  原始数据 → 过滤后数据 + token编码Arrow数据 + 细胞类型聚合数据
python preprocess_data.py \
    --h5ad ${RAW_H5AD} \
    --output data/${ORGAN} \
    --tokens prior_knowledge/human_mouse_tokens.pickle \
    --medians prior_knowledge/public/human_gene_median_after_filter.pickle
#  生成: data/${ORGAN}/filtered.h5ad              过滤后数据
#        data/${ORGAN}/single_cell_dataset/        token编码Arrow数据
#        data/${ORGAN}/cell_type_aggregated/       细胞类型聚合数据（Step 4用）

# ====== Step 1: CellChat 分析 ======
#  过滤后数据 → 细胞互作矩阵 + 通讯概率
python CellChatAnalysis/h5ad_to_csv.py \
    --input data/${ORGAN}/filtered.h5ad \
    --output data/${ORGAN}/cellchat/
Rscript CellChatAnalysis/csv_to_rds.R data/${ORGAN}/cellchat/ 3 200
Rscript CellChatAnalysis/cellchat_gold_standard.R data/${ORGAN}/cellchat/ 4
#  生成: data/${ORGAN}/cellchat/cell_interaction_strength_matrix.csv
#        data/${ORGAN}/cellchat/cellchat_communication.csv

# ====== Step 2: CellPhoneDB v5 分析 ======
#  原始数据 → 统计显著相互作用矩阵
CPDB_ZIP=/path/to/cellphonedb.zip
CPDB_DATA=/path/to/CellPhoneAnalysis/v5.0.0/
python run_cpdb.py \
    --h5ad ${RAW_H5AD} \
    --cpdb_db ${CPDB_ZIP} \
    --cpdb_genes ${CPDB_DATA} \
    --output data/${ORGAN}/cellphonedb/
#  生成: data/${ORGAN}/cellphonedb/significant_means.txt

# ====== Step 3: 构建联合金标准 ======
#  权重自动从来源显著性学习 (w ∝ 平均显著LR配对数)
python genecompass_gold_standard.py \
    --cellchat data/${ORGAN}/cellchat \
    --cpdb data/${ORGAN}/cellphonedb \
    --output data/${ORGAN}/gold_standard
#  生成: data/${ORGAN}/gold_standard/complete_labeled_interactions.csv

# ====== Step 4: 5折交叉验证 ======
#  用金标准标签训练GeneCompass，评估Spearman秩相关系数
python pipeline_cv.py \
    --proj_root . \
    --gs_path data/${ORGAN}/gold_standard/complete_labeled_interactions.csv \
    --dataset data/${ORGAN}/cell_type_aggregated \
    --output results/${ORGAN}_cv \
    --organ ${ORGAN^} --epochs 30 --batch 1 --grad_accum 4
#  生成: results/${ORGAN}_cv/cv_summary.json + fold{1-5}/ + 可视化

# ====== Step 5: 独立推理（可选） ======
python pipeline_inference.py \
    --model results/${ORGAN}_cv/fold1/best_model \
    --test_set results/${ORGAN}_cv/fold1/data_splits/test \
    --token_dict prior_knowledge/human_mouse_tokens.pickle
```

## 评估指标

### 主要指标：Spearman 秩相关系数 ρ
衡量模型对细胞类型互作强度排序与金标准共识排序的一致性。

| ρ 范围 | 等级 |
|---------|------|
| ρ ≥ 0.7 | 优秀 |
| 0.5 ≤ ρ < 0.7 | 良好 |
| 0.3 ≤ ρ < 0.5 | 中等 |
| ρ < 0.3 | 弱 |

每折报告：ρ + Bootstrap 95%置信区间 + 置换检验p值。
最终结果：5折均值 ± 标准差。

### 辅助指标
- **Pearson r**：与金标准的线性相关性
- **R²**：可解释方差
- **RMSE**：均方根误差

## 输出结构

```
data/{organ}/
├── filtered.h5ad                       # 过滤后数据 (Step 0)
├── single_cell_dataset/                # token编码Arrow数据 (Step 0)
├── cell_type_aggregated/               # 细胞类型聚合数据 (Step 0)
├── cellchat/                           # CellChat输出 (Step 1)
│   ├── cell_interaction_strength_matrix.csv
│   └── cellchat_communication.csv
├── cellphonedb/                        # CellPhoneDB输出 (Step 2)
│   ├── significant_means.txt
│   └── statistical_analysis_*.txt
└── gold_standard/                      # 联合金标准 (Step 3)
    ├── complete_labeled_interactions.csv
    └── gold_standard_stats.json

results/{organ}_cv/                     # 5折交叉验证输出 (Step 4)
├── cv_summary.json                     # 5折汇总：均值±标准差
├── fold{1-5}/
│   ├── metrics.json                    # 每折 ρ + CI + Pearson + R²
│   ├── best_model/                     # 训练模型 (pytorch_model.bin)
│   ├── test_true.npy                   # 真实标签
│   └── test_pred.npy                   # 预测值
├── interaction_heatmap.png             # 300dpi 可视化
├── interaction_network.png
├── interaction_circular.png
├── interaction_bubble.png
├── interaction_flow.png
├── autocrine_scores.png
└── true_vs_predicted.png
```

## 引用

```bibtex
@software{ccc_genecompass,
  title = {CCC-GeneCompass: Cell-Cell Communication via Large Language Model},
  year = {2025},
  note = {Based on GeneCompass: A Large-Scale Pretrained Model for Single-Cell Gene Expression}
}
```

## 许可

MIT
>>>>>>> 0a00b04 (Initial commit: CCC-GeneCompass v3 - gene compass BERT model for cell-cell communication analysis)
