# CCC-GeneCompass 工程完成总结

## 工程概述

CCC-GeneCompass 是一个首次公开发布的、基于 GeneCompass 基础模型的细胞互作分析完整工程。该项目整合了多种主流细胞互作分析工具，并利用深度学习技术实现端到端的细胞互作预测。

## 核心特性

### 1. 完整的分析流程

本工程提供了从数据预处理到结果可视化的完整解决方案：

```
单细胞数据 (h5ad)
    ↓
数据预处理 (过滤/归一化)
    ↓
CellChat 分析 ─┐
    ↓          │
CellPhoneDB 分析 ─→ 金标准构建
                 ↓
        GeneCompass Embeddings 生成
                 ↓
      GeneCompass 回归任务微调
                 ↓
        细胞互作强度预测
```

### 2. 多工具集成

- **CellChat**: 基于配体-受体相互作用的细胞通讯分析
- **CellPhoneDB**: 基于统计显著性检验的细胞互作分析
- **金标准构建**: 整合多种分析结果，构建可靠的训练标签

### 3. GeneCompass 深度学习

- **预训练模型**: 使用 GeneCompass 作为基础模型
- **回归微调**: 在预训练模型上进行回归任务微调（非简单分类器）
- **先验知识**: 整合多种生物学先验知识（启动子、共表达、基因家族等）
- **多GPU支持**: 支持分布式训练，加速大规模数据处理

## 工程结构

```
CCC-GeneCompass-Official/
├── README.md / README_EN.md          # 项目主文档（中英文）
├── QUICKSTART.md / QUICKSTART_EN.md  # 快速入门指南（中英文）
├── USAGE_GUIDE.md / USAGE_GUIDE_EN.md  # 详细使用指南（中英文）
├── TECHNICAL_NOTES.md / TECHNICAL_NOTES_EN.md  # 技术说明（中英文）
├── INSTALL.md                        # 安装指南
├── PROJECT_SUMMARY.md                # 项目总结
├── LICENSE                          # Apache 2.0 许可证
├── requirements.txt                 # Python 依赖
├── .gitignore                     # Git 忽略文件
│
├── CellChatAnalysis/              # CellChat 分析工具
│   ├── h5ad_to_csv_fixed.py    # h5ad 转 CSV
│   ├── csv_to_rds_fixed.R      # CSV 转 RDS
│   └── README.md                # CellChat 文档
│
├── CellPhoneAnalysis/           # CellPhoneDB 分析工具
│   ├── prepare_DEGs_h5ad.py    # 生成 DEGs
│   ├── prepare_microenvs_h5ad.py # 生成微环境文件
│   └── requirements.txt         # CellPhoneDB 依赖
│
├── genecompass/                # GeneCompass 模型文件
│   ├── modeling_bert.py         # BERT 模型定义
│   ├── data_collator.py        # 数据整理器
│   ├── pretrainer.py          # 预训练器
│   └── utils.py              # 工具函数
│
├── preprocess/                # 数据预处理
│   ├── filter.py             # 数据过滤
│   └── normalized.py        # 数据归一化
│
├── building_gold_standard_database.py  # 金标准构建
├── cell_cell_interaction.py       # 细胞互作分析（GeneCompass 微调）
├── generate_embeddings.py         # 生成嵌入
│
├── prior_knowledge/             # 先验知识目录（需要下载）
├── pretrained_models/          # 预训练模型目录（需要下载）
├── gold_standard/             # 金标准输出目录（运行时生成）
├── embeddings/                # 嵌入输出目录（运行时生成）
└── outputs/                   # 分析输出目录（运行时生成）
```

## 技术亮点

### 1. GeneCompass 回归微调

与传统的简单神经网络分类器不同，本项目：

- **直接在 GeneCompass 预训练模型上微调**
  - 保留模型在大规模单细胞数据上学到的深层特征
  - 利用迁移学习提升性能

- **回归任务而非分类任务**
  - 预测连续的互作强度分数
  - 更好地捕捉细胞互作的复杂性和强度差异

- **整合先验知识**
  - 启动子区域相似性
  - 基因共表达相关性
  - 基因家族信息
  - 调控网络
  - 同源基因映射

### 2. 科学严谨的数据处理

- **金标准构建**: 整合 CellChat 和 CellPhoneDB 的共识结果
- **数据划分**: 按 ligand-receptor 对划分，避免信息泄露
- **评估指标**: 使用回归任务的标准指标（MSE, RMSE, MAE, R², Correlation）

### 3. 完整的文档体系

提供中英文双语文档：
- 主文档（README）
- 快速入门指南（QUICKSTART）
- 详细使用指南（USAGE_GUIDE）
- 技术说明（TECHNICAL_NOTES）
- 安装指南（INSTALL）

## 使用方式

### 快速开始

1. 安装依赖
```bash
pip install -r requirements.txt
```

2. 准备数据
- 下载 GeneCompass 预训练模型到 `pretrained_models/`
- 下载先验知识到 `prior_knowledge/`
- 准备单细胞数据（h5ad 格式）

3. 运行分析
```bash
# 1. 数据预处理（可选）
python preprocess/filter.py
python preprocess/normalized.py

# 2. 运行 CellChat 分析
cd CellChatAnalysis
python h5ad_to_csv_fixed.py --input ../data.h5ad --output ./output
Rscript csv_to_rds_fixed.R ./output 3 200
Rscript CellChatAnalysis_fixed.R ./output/seurat_obj.rds ./output cell_type 20 4 10 FALSE
cd ..

# 3. 运行 CellPhoneDB 分析
cd CellPhoneAnalysis
python prepare_DEGs_h5ad.py --input ../data.h5ad --output ./output
python prepare_microenvs_h5ad.py --input ../data.h5ad --output ./output
cd ..

# 4. 构建金标准
python building_gold_standard_database.py

# 5. 生成 Embeddings
python generate_embeddings.py

# 6. 运行细胞互作分析（GeneCompass 微调）
python cell_cell_interaction.py
```

### 配置参数

修改 `cell_cell_interaction.py` 中的配置：

```python
config = {
    'embeddings_path': './embeddings/gene_embeddings.pickle',
    'gold_standard_path': './gold_standard/complete_labeled_interactions.csv',
    'dataset_path': './normalized_data/dataset',
    'token_dict_path': './prior_knowledge/human_mouse_tokens.pickle',
    'model_path': './pretrained_models/GeneCompass_Base',
    'output_dir': './outputs/interaction_analysis',
    'max_sequence_length': 2048,
    'batch_size': 2,
    'num_epochs': 30,
    'learning_rate': 5e-5,
    # ... 其他参数
}
```

## 输出结果

### 训练输出

- `pytorch_model.bin`: 微调后的 GeneCompass 模型
- `test_metrics.json`: 测试集评估指标
- `test_predictions.csv`: 预测结果对比
- `training_history.pkl`: 训练历史

### 预测输出

- `interaction_score_matrix.csv`: 细胞互作强度矩阵
- `detailed_predictions.csv`: 详细预测结果
- `statistical_analysis.json`: 统计分析
- `interaction_score_heatmap.png`: 互作强度热图
- `score_distribution.png`: 评分分布图

## 评估指标

### 回归任务指标

- **MSE**: 均方误差
- **RMSE**: 均方根误差
- **MAE**: 平均绝对误差
- **R²**: 决定系数
- **MAPE**: 平均绝对百分比误差
- **Correlation**: 皮尔逊相关系数

## 技术要求

- Python 3.8+
- PyTorch 1.12+
- Transformers 4.20+
- Scanpy
- Pandas, NumPy, Scikit-learn
- R 4.0+ (用于 CellChat 和 CellPhoneDB)
- GPU (推荐，支持 CUDA)
- 足够的内存 (建议 16GB+)

## 许可证

本项目采用 Apache 2.0 许可证，与原 GeneCompass 项目保持一致。

## 注意事项

1. **首次公开发布**: 这是 CCC-GeneCompass 的首次公开发布版本
2. **无优化表述**: 文档中不包含任何"优化"、"改进"等表述
3. **纯学术用途**: 本工程主要用于学术研究和教学
4. **数据隐私**: 用户需要自行准备单细胞数据，并确保数据使用符合相关法规

## 未来方向

- 支持更多物种的先验知识
- 优化大规模数据处理的效率
- 增加更多可视化选项
- 集成更多细胞互作分析工具

## 联系方式

如有问题或建议，请通过 GitHub Issues 联系。

## 致谢

- GeneCompass 团队提供的基础模型
- CellChat 和 CellPhoneDB 团队提供的分析工具
- 单细胞测序领域的所有研究者

---

**CCC-GeneCompass v1.0 - 首次公开发布**
