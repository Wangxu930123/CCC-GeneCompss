# 技术说明

## 核心技术特点

### 基于GeneCompass的回归任务微调

本项目采用基于GeneCompass基础模型的回归任务微调方法进行细胞互作分析，而非简单的分类模型。

#### 关键技术优势

1. **直接微调预训练模型**
   - 在GeneCompass预训练模型基础上进行微调
   - 保留模型在大规模单细胞数据上学到的深层特征表示
   - 避免从零训练小模型带来的特征表达能力损失

2. **回归任务设计**
   - 预测连续的细胞互作强度分数
   - 使用金标准中的共识分数（Consensus_Score）作为训练标签
   - 更好地捕捉细胞互作的复杂性和强度差异

3. **多GPU分布式训练**
   - 支持单机多GPU训练
   - 使用PyTorch分布式训练框架
   - 自动处理数据同步和模型并行

#### 模型架构

```
GeneCompass预训练模型
├── BERT编码器（多层Transformer）
│   ├── Self-attention机制
│   ├── 前馈神经网络
│   └── 残差连接和层归一化
├── 先验知识注入
│   ├── 启动子区域相似性
│   ├── 共表达相关性
│   ├── 基因家族信息
│   ├── 调控网络
│   └── 同源基因映射
└── 回归预测头
    └── 单输出节点（互作强度分数）
```

#### 训练流程

1. **数据准备**
   - 从金标准中提取细胞对和共识分数
   - 构建细胞对序列（发送者序列 + 接收者序列）
   - 加载预训练的GeneCompass模型和先验知识

2. **模型微调**
   - 使用HuggingFace Transformers框架
   - 回归损失函数：MSELoss
   - 学习率调度：带预热的线性衰减
   - 评估指标：MSE, RMSE, MAE, R², Correlation

3. **预测与可视化**
   - 使用微调后的模型预测细胞互作分数
   - 生成互作强度矩阵
   - 绘制热图和分布图

#### 与传统方法的区别

| 特性 | 传统分类方法 | GeneCompass回归微调 |
|------|-------------|-------------------|
| 模型基础 | 简单神经网络 | 预训练大模型 |
| 任务类型 | 二分类 | 回归 |
| 输出 | 互作概率 | 互作强度分数 |
| 特征学习 | 从零开始 | 迁移学习 |
| 先验知识 | 无 | 整合多种生物学先验 |
| 预测粒度 | 有/无互作 | 连续强度值 |

#### 性能优势

- **更好的特征表示**：继承GeneCompass在大规模数据上学到的基因表示
- **更准确的预测**：回归任务能更好地建模互作强度的连续性
- **更强的泛化能力**：预训练模型提供良好的初始化
- **生物学可解释性**：整合多种先验知识，增强可解释性

## 适用场景

- 单细胞RNA测序数据的细胞互作分析
- 需要精确预测互作强度的研究
- 跨物种细胞互作比较
- 大规模单细胞数据集的高效分析

## 技术要求

- Python 3.8+
- PyTorch 1.12+
- Transformers 4.20+
- GPU（推荐，支持CUDA）
- 足够的内存（建议16GB+）

## 参考文献

- GeneCompass: A foundation model for single-cell RNA sequencing analysis
- CellChat: Inferring and analyzing cell-cell communication
- CellPhoneDB: Inferring cell-cell communication from combined expression of multi-subunit ligand-receptor complexes
