# CCC-GeneCompass Project Summary

## 项目概述 / Project Overview

CCC-GeneCompass 是一个基于 GeneCompass 单细胞基础模型的细胞互作分析工具。本项目首次公开，整合了多种细胞通信分析方法，提供了从数据预处理到互作预测的完整流程。

CCC-GeneCompass is a cell-cell interaction analysis tool based on the GeneCompass single-cell foundation model. This project is being released for the first time, integrating multiple cell communication analysis methods to provide a complete pipeline from data preprocessing to interaction prediction.

## 核心功能 / Core Features

### 1. 细胞互作分析工具集成 / Cell-Cell Interaction Analysis Tool Integration

- **CellChat**: 基于配体-受体相互作用的概率模型
- **CellPhoneDB**: 基于统计显著性检验的分析方法

### 2. 金标准构建 / Gold Standard Construction

- 整合 CellChat 和 CellPhoneDB 的分析结果
- 使用共识分数评估互作置信度
- 支持多种阈值方法（分位数/绝对阈值）
- 生成完整的机器学习数据集

### 3. GeneCompass 模型应用 / GeneCompass Model Application

- 基于预训练的 GeneCompass 模型生成细胞/基因嵌入
- 使用金标准微调模型进行细胞互作预测
- 支持端到端的单细胞数据到互作网络分析

## 技术架构 / Technical Architecture

### 分析流程 / Analysis Pipeline

```
单细胞数据 (.h5ad)
    ↓
数据预处理 (过滤、归一化)
    ↓
    ├─→ CellChat 分析
    │      - 转换 h5ad 为 Seurat 对象
    │      - 计算互作概率
    │      - 生成互作矩阵
    │
    ├─→ CellPhoneDB 分析
    │      - 生成差异表达基因
    │      - 统计显著性检验
    │      - 计算互作强度
    │
    ↓
金标准构建
    - 归一化分数
    - 计算共识分数
    - 分配标签
    ↓
生成 Embeddings
    - 加载 GeneCompass 模型
    - 提取细胞嵌入
    ↓
细胞互作预测
    - 训练神经网络
    - 预测互作
    ↓
结果可视化
```

### 核心算法 / Core Algorithms

#### 共识分数计算 / Consensus Score Calculation

```
Consensus_Score = (Norm_CellChat + Norm_CPDB_Mean + Norm_CPDB_Max) / 3
```

其中：
- `Norm_CellChat`: 归一化后的 CellChat 互作分数
- `Norm_CPDB_Mean`: 归一化后的 CellPhoneDB 平均互作分数
- `Norm_CPDB_Max`: 归一化后的 CellPhoneDB 最大互作分数

#### 归一化方法 / Normalization Method

使用 MinMaxScaler 将所有分数归一化到 [0, 1] 范围。

## 文件结构 / File Structure

```
CCC-GeneCompass-Official/
├── README.md                       # 主要文档（中文）
├── README_EN.md                    # 主要文档（英文）
├── QUICKSTART.md                   # 快速入门指南（中文）
├── QUICKSTART_EN.md                # 快速入门指南（英文）
├── INSTALL.md                      # 安装指南（中英文）
├── LICENSE                         # Apache 2.0 许可证
├── requirements.txt                # Python 依赖
├── .gitignore                     # Git 忽略文件
│
├── CellChatAnalysis/              # CellChat 分析工具
│   ├── h5ad_to_csv_fixed.py      # h5ad 转 CSV
│   ├── csv_to_rds_fixed.R        # CSV 转 RDS
│   ├── CellChatAnalysis_fixed.R    # CellChat 主脚本
│   └── README.md                  # CellChat 文档
│
├── CellPhoneAnalysis/             # CellPhoneDB 分析工具
│   ├── prepare_DEGs_h5ad.py      # 生成 DEGs
│   ├── prepare_microenvs_h5ad.py # 生成微环境文件
│   ├── CellPhoneAnalysis.py       # CellPhoneDB 主脚本
│   ├── requirements.txt           # CellPhoneDB 依赖
│   └── ANALYSIS_REPORT.md       # 分析报告
│
├── genecompass/                   # GeneCompass 模型文件
│   ├── modeling_bert.py          # BERT 模型定义
│   ├── data_collator.py         # 数据整理器
│   ├── pretrainer.py            # 预训练器
│   └── utils.py                # 工具函数
│
├── preprocess/                    # 数据预处理
│   ├── filter.py                # 数据过滤
│   └── normalized.py           # 数据归一化
│
├── building_gold_standard_database.py  # 金标准构建
├── cell_cell_interaction.py       # 细胞互作分析
├── generate_embeddings.py         # 生成嵌入
│
├── gold_standard/                # 金标准输出目录（创建时生成）
├── embeddings/                   # 嵌入输出目录（创建时生成）
├── outputs/                      # 分析输出目录（创建时生成）
└── pretrained_models/             # 预训练模型目录（需要手动下载）
```

## 依赖关系 / Dependencies

### Python 依赖 / Python Dependencies

- `torch>=1.13.1`: 深度学习框架
- `transformers==4.30.0`: Hugging Face Transformers
- `scanpy>=1.9`: 单细胞数据分析
- `anndata>=0.9`: AnnData 数据结构
- `pandas>=2.0`: 数据处理
- `numpy>=1.23`: 数值计算
- `scikit-learn>=1.2`: 机器学习
- `matplotlib>=3.7`: 可视化
- `seaborn>=0.12`: 统计可视化
- `networkx>=2.8`: 网络分析

### R 依赖 / R Dependencies

- `Seurat`: 单细胞数据分析
- `CellChat`: 细胞通信分析
- `dplyr`: 数据处理
- `ggplot2`: 可视化
- `ComplexHeatmap`: 热图
- `patchwork`: 图形组合

## 使用示例 / Usage Examples

### 完整流程示例 / Complete Pipeline Example

```bash
# 1. CellChat 分析
cd CellChatAnalysis
python h5ad_to_csv_fixed.py --input ../data.h5ad --output ./output
Rscript csv_to_rds_fixed.R ./output
Rscript CellChatAnalysis_fixed.R ./output/seurat_obj.rds ./output

# 2. CellPhoneDB 分析
cd ../CellPhoneAnalysis
python prepare_DEGs_h5ad.py --h5ad ../data.h5ad --outdir ./output
python prepare_microenvs_h5ad.py --h5ad ../data.h5ad --output ./output/microenv.tsv
python CellPhoneAnalysis.py --h5ad ../data.h5ad --outdir ./output

# 3. 构建金标准
cd ..
python building_gold_standard_database.py \
  --cellchat_dir ./CellChatAnalysis/output/ \
  --cpdb_dir ./CellPhoneAnalysis/output/ \
  --output_dir ./gold_standard

# 4. 生成嵌入
python generate_embeddings.py \
  --dataset_path ./normalized_data/ \
  --model_path ./pretrained_models/GeneCompass_Base \
  --token_dict_path ./prior_knowledge/human_mouse_tokens.pickle \
  --output_path ./embeddings/embeddings.pickle

# 5. 细胞互作分析
python cell_cell_interaction.py
```

## 输出说明 / Output Description

### CellChat 输出 / CellChat Outputs

- `cellchat_communication.csv`: 细胞通讯事件
- `cellchat_pathways.csv`: 信号通路
- `cell_interaction_strength_matrix.csv`: 互作强度矩阵
- `communication_network.png`: 通讯网络图
- `communication_heatmap.png`: 互作热图

### CellPhoneDB 输出 / CellPhoneDB Outputs

- `significant_means.txt`: 显著互作结果
- `means.txt`: 所有互作均值
- `pvalues.txt`: 显著性检验结果

### 金标准输出 / Gold Standard Outputs

- `complete_labeled_interactions.csv`: 完整标记数据集
- `machine_learning_dataset.csv`: 机器学习数据集
- `gold_standard_interactions.csv`: 金标准正样本
- `dataset_statistics.csv`: 数据集统计
- 多种可视化图表

### GeneCompass 输出 / GeneCompass Outputs

- `embeddings.pickle`: 细胞/基因嵌入向量
- `best_model.pt`: 训练好的模型
- `evaluation_results.pkl`: 评估结果

## 项目特点 / Project Highlights

1. **首次公开 / First Public Release**: 这是 CCC-GeneCompass 的首次公开发布
2. **完整流程 / Complete Pipeline**: 从数据预处理到结果可视化的完整解决方案
3. **多工具整合 / Multi-tool Integration**: 整合 CellChat 和 CellPhoneDB 两种主流工具
4. **深度学习 / Deep Learning**: 利用 GeneCompass 基础模型进行预测
5. **模块化设计 / Modular Design**: 每个模块可独立使用
6. **详细文档 / Detailed Documentation**: 提供中英文双语文档

## 未来方向 / Future Directions

- 支持更多细胞互作分析工具
- 优化模型性能和准确性
- 增加更多可视化选项
- 支持分布式训练
- 提供预训练的分析模型

## 引用 / Citation

如果您在研究中使用了本工具，请引用 / If you use this tool in your research, please cite:

- **GeneCompass**: [论文链接] / [Paper Link]
- **CellChat**: Jin, S., et al. (2021). Nature Protocols.
- **CellPhoneDB**: Vento-Tormo, R., et al. (2018). Science.

## 许可证 / License

本项目遵循 Apache 2.0 许可证 / This project follows the Apache 2.0 License.

## 联系方式 / Contact

如有问题或建议，请提交 Issue 或 Pull Request。
For questions or suggestions, please submit an Issue or Pull Request.

---

**版本 / Version**: 1.0
**发布日期 / Release Date**: March 2026
**状态 / Status**: 首次公开发布 / First Public Release
