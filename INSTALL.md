# Installation Guide / 安装指南

## English / 英文

### System Requirements / 系统要求

- **Operating System**: Linux, macOS, or Windows
- **Python**: 3.8 or higher
- **R**: 4.0 or higher
- **Memory**: At least 16GB RAM (recommended 32GB+)
- **Storage**: At least 20GB free space
- **GPU**: Optional (recommended for faster training)

### Python Environment Setup / Python 环境设置

```bash
# Create conda environment / 创建 conda 环境
conda create -n ccc python=3.10
conda activate ccc

# Install Python dependencies / 安装 Python 依赖
pip install -r requirements.txt
```

If you encounter issues with `transformers==4.30.0`:

```bash
# Install tokenizers first / 先安装 tokenizers
conda install -c conda-forge tokenizers=0.13.3
pip install transformers==4.30.0
```

### R Environment Setup / R 环境设置

```r
# Install CRAN packages / 安装 CRAN 包
install.packages(c("Seurat", "ggplot2", "patchwork", "dplyr",
                   "future", "ComplexHeatmap", "RColorBrewer", "stringr"))

# Install Bioconductor packages / 安装 Bioconductor 包
if (!requireNamespace("BiocManager", quietly = TRUE))
    install.packages("BiocManager")
BiocManager::install("ComplexHeatmap")

# Install CellChat / 安装 CellChat
install.packages("devtools")
devtools::install_github("sqjin/CellChat")
```

Alternatively, use conda:

```bash
conda install -c bioconda r-seurat r-cellchat
```

### CellPhoneDB Installation / CellPhoneDB 安装

For detailed CellPhoneDB installation, please refer to the official documentation:
https://github.com/ventolab/CellPhoneDB

## 中文

### 系统要求

- **操作系统**: Linux、macOS 或 Windows
- **Python**: 3.8 或更高版本
- **R**: 4.0 或更高版本
- **内存**: 至少 16GB RAM（推荐 32GB+）
- **存储**: 至少 20GB 可用空间
- **GPU**: 可选（推荐用于加速训练）

### Python 环境设置

```bash
# 创建 conda 环境
conda create -n ccc python=3.10
conda activate ccc

# 安装 Python 依赖
pip install -r requirements.txt
```

如果 `transformers==4.30.0` 安装遇到问题：

```bash
# 先安装 tokenizers
conda install -c conda-forge tokenizers=0.13.3
pip install transformers==4.30.0
```

### R 环境设置

```r
# 安装 CRAN 包
install.packages(c("Seurat", "ggplot2", "patchwork", "dplyr",
                   "future", "ComplexHeatmap", "RColorBrewer", "stringr"))

# 安装 Bioconductor 包
if (!requireNamespace("BiocManager", quietly = TRUE))
    install.packages("BiocManager")
BiocManager::install("ComplexHeatmap")

# 安装 CellChat
install.packages("devtools")
devtools::install_github("sqjin/CellChat")
```

或者使用 conda：

```bash
conda install -c bioconda r-seurat r-cellchat
```

### CellPhoneDB 安装

详细的 CellPhoneDB 安装请参考官方文档：
https://github.com/ventolab/CellPhoneDB

## Data Preparation / 数据准备

### Download Pretrained Models / 下载预训练模型

1. Visit GeneCompass official repository: https://github.com/xCompass-AI/GeneCompass
2. Download pretrained models:
   - GeneCompass_Base (recommended for production)
   - GeneCompass_Small (for testing)
3. Place models in `./pretrained_models/` directory

### Download Prior Knowledge / 下载先验知识

Download the following files from GeneCompass repository and place them in `./prior_knowledge/`:

- `human_mouse_tokens.pickle`
- `human_gene_median_after_filter.pickle`
- Other required knowledge files

### Prepare Your Data / 准备您的数据

Ensure your single-cell data is in `.h5ad` format and contains:
- Expression matrix (`adata.X`)
- Cell metadata (`adata.obs`) with cell type column
- Gene information (`adata.var`)

## Verification / 验证

After installation, verify your setup:

```bash
# Test Python environment / 测试 Python 环境
python -c "import torch; import transformers; import scanpy; print('Python environment OK')"

# Test R environment (in R) / 测试 R 环境
R
> library(Seurat)
> library(CellChat)
> print('R environment OK')
```

## Troubleshooting / 故障排除

### Issue: Module not found / 问题: 模块未找到

```bash
# Reinstall the missing package / 重新安装缺失的包
pip install <package_name>
```

### Issue: R package installation fails / 问题: R 包安装失败

```r
# Try installing from a different source / 尝试从不同源安装
install.packages("<package_name>", repos="https://cran.r-project.org")
```

### Issue: Out of memory / 问题: 内存不足

- Reduce batch size in configuration / 减小配置中的批次大小
- Use a smaller dataset for testing / 使用较小的数据集进行测试
- Close other applications / 关闭其他应用程序

## Next Steps / 下一步

After successful installation, refer to:
- [QUICKSTART.md](QUICKSTART.md) - Quick start guide (Chinese) / 快速入门指南（中文）
- [QUICKSTART_EN.md](QUICKSTART_EN.md) - Quick start guide (English) / 快速入门指南（英文）
- [README.md](README.md) - Main documentation (Chinese) / 主要文档（中文）
- [README_EN.md](README_EN.md) - Main documentation (English) / 主要文档（英文）
