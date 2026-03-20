# CCC-GeneCompass 使用指南 / Usage Guide

## 首次使用快速检查 / First-Time Setup Checklist

### 1. 环境检查 / Environment Check

```bash
# Python 环境 / Python Environment
python --version  # Should be 3.8+
conda info --envs

# 检查依赖 / Check dependencies
python -c "import torch; import transformers; import scanpy; print('✓ Python OK')"

# R 环境 / R Environment
R --version  # Should be 4.0+

# 在 R 中检查 / In R, check:
library(Seurat)
library(CellChat)
print('✓ R OK')
```

### 2. 数据准备 / Data Preparation

**必需文件 / Required Files:**

- [ ] 预训练模型在 `./pretrained_models/` 目录下
- [ ] Token 字典在 `./prior_knowledge/` 目录下
- [ ] 单细胞数据 (`.h5ad` 格式）准备好
- [ ] 数据包含细胞类型信息 (`cell_type` 列）

### 3. 目录结构验证 / Directory Structure Verification

```bash
# 运行此命令检查目录结构 / Run this to check directory structure
ls -la
ls CellChatAnalysis/
ls CellPhoneAnalysis/
ls genecompass/
ls preprocess/
```

## 完整分析流程 / Complete Analysis Workflow

### 阶段 1: 数据预处理 / Phase 1: Data Preprocessing

**可选步骤 / Optional Step** - 如果您的数据已经归一化，可以跳过

```bash
# 步骤 1.1: 数据过滤 / Step 1.1: Data Filtering
cd preprocess
# 修改 filter.py 中的路径指向您的数据
python filter.py

# 步骤 1.2: 数据归一化 / Step 1.2: Data Normalization
# 修改 normalized.py 中的路径指向您的数据
python normalized.py
cd ..
```

**预期输出 / Expected Output:**
- `./filtered_data/` - 过滤后的数据
- `./normalized_data/` - 归一化后的数据（包含 dataset 文件夹）

### 阶段 2: CellChat 分析 / Phase 2: CellChat Analysis

```bash
cd CellChatAnalysis

# 步骤 2.1: h5ad 转 CSV / Step 2.1: Convert h5ad to CSV
python h5ad_to_csv_fixed.py \
  --input ../path/to/your_data.h5ad \
  --output ./cellchat_output \
  --celltype_col cell_type

# 步骤 2.2: CSV 转 Seurat 对象 / Step 2.2: Convert CSV to Seurat
Rscript csv_to_rds_fixed.R ./cellchat_output 3 200

# 步骤 2.3: 运行 CellChat 分析 / Step 2.3: Run CellChat
Rscript CellChatAnalysis_fixed.R \
  ./cellchat_output/seurat_obj.rds \
  ./cellchat_output \
  cell_type \
  20 \
  4 \
  10 \
  FALSE

cd ..
```

**预期输出 / Expected Output:**
- `./CellChatAnalysis/cellchat_output/cellchat_communication.csv`
- `./CellChatAnalysis/cellchat_output/cellchat_pathways.csv`
- `./CellChatAnalysis/cellchat_output/cell_interaction_strength_matrix.csv`
- 多种可视化图表

**常见问题 / Common Issues:**
- 如果 R 脚本报错 "multicore not supported"，脚本会自动检测并使用 "multisession"
- 如果基因名是 Ensembl ID，需要先转换为 Gene Symbol
- 如果没有检测到互作，检查细胞类型数量是否足够（至少 2 种）

### 阶段 3: CellPhoneDB 分析 / Phase 3: CellPhoneDB Analysis

```bash
cd CellPhoneAnalysis

# 步骤 3.1: 生成 DEGs / Step 3.1: Generate DEGs
python prepare_DEGs_h5ad.py \
  --mode degs_only \
  --h5ad ../path/to/your_data.h5ad \
  --outdir ./cpdb_output \
  --groupby cell_type \
  --n_top_genes 250

# 步骤 3.2: 生成微环境文件 / Step 3.2: Generate microenv file
python prepare_microenvs_h5ad.py \
  --h5ad ../path/to/your_data.h5ad \
  --output ./cpdb_output/microenv.tsv \
  --groupby cell_type

# 步骤 3.3: 运行 CellPhoneDB / Step 3.3: Run CellPhoneDB
python CellPhoneAnalysis.py \
  --h5ad ../path/to/your_data.h5ad \
  --cpdb ./path/to/cellphonedb.zip \
  --degs ./cpdb_output/DEGs.tsv \
  --microenv ./cpdb_output/microenv.tsv \
  --groupby cell_type \
  --outdir ./cpdb_output

cd ..
```

**预期输出 / Expected Output:**
- `./CellPhoneAnalysis/cpdb_output/significant_means.txt`
- `./CellPhoneAnalysis/cpdb_output/means.txt`
- `./CellPhoneAnalysis/cpdb_output/pvalues.txt`

**注意 / Note:**
- 需要先安装 CellPhoneDB，参考 CellPhoneDB 官方文档
- 确保 CellPhoneDB 数据库文件路径正确

### 阶段 4: 构建金标准 / Phase 4: Build Gold Standard

```bash
# 步骤 4: 构建金标准数据集 / Step 4: Build gold standard dataset
python building_gold_standard_database.py \
  --cellchat_dir ./CellChatAnalysis/cellchat_output/ \
  --cpdb_dir ./CellPhoneAnalysis/cpdb_output/ \
  --output_dir ./gold_standard/ \
  --threshold_method quantile \
  --threshold_value 0.7
```

**预期输出 / Expected Output:**
- `./gold_standard/complete_labeled_interactions.csv`
- `./gold_standard/machine_learning_dataset.csv`
- `./gold_standard/gold_standard_interactions.csv`
- `./gold_standard/dataset_statistics.csv`
- 多种可视化图表

**参数说明 / Parameter Explanation:**
- `--threshold_method quantile`: 使用分位数阈值（推荐）
- `--threshold_value 0.7`: 前面 30% 的互作作为金标准正样本
- 可以调整为其他值（如 0.5 或 0.8）来改变正样本数量

### 阶段 5: 生成 Embeddings / Phase 5: Generate Embeddings

```bash
# 步骤 5: 生成 GeneCompass 嵌入 / Step 5: Generate GeneCompass embeddings
python generate_embeddings.py \
  --dataset_path ./normalized_data/TabulaSapiens/tabula_sapiens_liver/ \
  --model_path ./pretrained_models/GeneCompass_Base \
  --token_dict_path ./prior_knowledge/human_mouse_tokens.pickle \
  --output_path ./embeddings/embeddings.pickle \
  --batch_size 128 \
  --gpu_ids "0"
```

**预期输出 / Expected Output:**
- `./embeddings/embeddings.pickle` - 嵌入向量文件

**注意 / Note:**
- 此步骤可能需要较长时间（数小时），取决于数据集大小
- 如果内存不足，减小 `--batch_size`
- 如果没有 GPU，去掉 `--gpu_ids` 参数或设置为空字符串

### 阶段 6: 细胞互作分析 / Phase 6: Cell-Cell Interaction Analysis

```bash
# 步骤 6: 训练和评估模型 / Step 6: Train and evaluate model
# 需要修改 cell_cell_interaction.py 中的路径
python cell_cell_interaction.py
```

**预期输出 / Expected Output:**
- `./outputs/best_model.pt` - 最佳模型
- `./outputs/evaluation_results.pkl` - 评估结果
- 训练日志显示每个 epoch 的性能

**配置参数 / Configuration:**
在 `cell_cell_interaction.py` 中修改以下变量：
```python
embeddings_path = './embeddings/embeddings.pickle'
gold_standard_path = './gold_standard/machine_learning_dataset.csv'
dataset_path = './normalized_data/TabulaSapiens/tabula_sapiens_liver/'
token_dict_path = './prior_knowledge/human_mouse_tokens.pickle'
output_dir = './outputs/'
```

## 结果解释 / Result Interpretation

### 金标准数据集 / Gold Standard Dataset

**文件结构 / File Structure:**
```csv
Sender,Receiver,Pair_ID,Gold_Standard_Label,Confidence_Level,Consensus_Score,...
CellA,CellB,CellA_CellB,1,Gold Standard,0.85,...
```

**列说明 / Column Description:**
- `Sender`: 发送细胞类型
- `Receiver`: 接收细胞类型
- `Gold_Standard_Label`: 1=金标准正样本，0=负样本
- `Confidence_Level`: 置信度级别
- `Consensus_Score`: 共识分数 (0-1)

### 评估指标 / Evaluation Metrics

**回归任务指标 / Regression Task Metrics:**
- **MSE (Mean Squared Error)**: 均方误差
- **RMSE (Root Mean Squared Error)**: 均方根误差
- **MAE (Mean Absolute Error)**: 平均绝对误差
- **R² (R-Squared)**: 决定系数
- **MAPE (Mean Absolute Percentage Error)**: 平均绝对百分比误差
- **Correlation**: 皮尔逊相关系数

**关键说明 / Key Notes:**
- 本项目使用基于GeneCompass的回归任务微调，而非简单的分类任务
- 模型直接在GeneCompass预训练模型上进行微调，保留其深层特征学习能力
- 预测输出为连续的互作强度分数，而非二元标签

## 常见问题排查 / Common Issues Troubleshooting

### 问题 1: 内存不足 / Issue 1: Out of Memory

**症状 / Symptoms:**
- Python: `MemoryError` 或程序崩溃
- R: 内存不足错误

**解决方案 / Solutions:**
```bash
# Python: 减小批次大小 / Python: Reduce batch size
python generate_embeddings.py --batch_size 64  # 从 128 减小到 64

# R: 减小并行 workers / R: Reduce parallel workers
Rscript CellChatAnalysis_fixed.R ... 2  # 从 4 减小到 2

# 关闭其他程序 / Close other applications
```

### 问题 2: GPU 未使用 / Issue 2: GPU Not Used

**症状 / Symptoms:**
- 训练速度很慢
- `nvidia-smi` 显示 GPU 利用率为 0%

**解决方案 / Solutions:**
```bash
# 检查 CUDA 可用性 / Check CUDA availability
python -c "import torch; print(torch.cuda.is_available())"

# 检查 PyTorch 版本 / Check PyTorch version
python -c "import torch; print(torch.__version__)"

# 如果返回 False，需要安装支持 CUDA 的 PyTorch 版本
```

### 问题 3: 基因名不匹配 / Issue 3: Gene Name Mismatch

**症状 / Symptoms:**
- CellChat 或 CellPhoneDB 没有检测到任何互作
- "0 features" 或类似错误

**解决方案 / Solutions:**
```bash
# 检查基因名格式 / Check gene name format
python -c "import scanpy as sc; adata=sc.read_h5ad('data.h5ad'); print(adata.var_names[:10])"

# 如果是 Ensembl ID (如 ENSG00000000003)，需要转换
# 使用 CellChatAnalysis 目录下的工具
```

### 问题 4: 路径错误 / Issue 4: Path Errors

**症状 / Symptoms:**
- `FileNotFoundError` 或文件不存在错误

**解决方案 / Solutions:**
```bash
# 使用绝对路径 / Use absolute paths
python script.py --input "/full/path/to/file.h5ad"

# 或者在代码中设置工作目录 / Or set working directory in code
import os
os.chdir('/path/to/project')
```

## 进阶使用 / Advanced Usage

### 自定义阈值 / Custom Thresholds

```bash
python building_gold_standard_database.py \
  --threshold_method absolute \
  --threshold_value 0.5  # 使用绝对阈值
```

### 使用不同的模型 / Use Different Models

```bash
python generate_embeddings.py \
  --model_path ./pretrained_models/GeneCompass_Small  # 使用小模型
```

### 调整模型架构 / Adjust Model Architecture

修改 `cell_cell_interaction.py` 中的模型定义：

```python
class SimpleCellInteractionClassifier(nn.Module):
    def __init__(self, input_dim, hidden_dim=512, num_classes=2):  # 增加隐藏层维度
        ...
```

## 下一步 / Next Steps

- 尝试在您自己的数据集上运行完整流程
- 根据结果调整参数（阈值、模型架构等）
- 探索不同的可视化选项
- 与传统方法结果进行比较

## 获取帮助 / Getting Help

1. 查看 [README.md](README.md) 了解项目概述
2. 查看 [QUICKSTART.md](QUICKSTART.md) 快速入门
3. 查看 [INSTALL.md](INSTALL.md) 安装指南
4. 提交 Issue 获取技术支持

---

**版本 / Version**: 1.0
**最后更新 / Last Updated**: March 2026
