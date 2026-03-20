# CCC-GeneCompass Quick Start Guide

## Overview

CCC-GeneCompass is a cell-cell interaction analysis tool based on the GeneCompass single-cell foundation model. This guide will help you quickly get started with the tool for cell-cell interaction analysis.

## Analysis Pipeline Overview

```
Single-cell data (.h5ad)
    │
    ├─────────────┬─────────────┬─────────────┐
    ▼             ▼             ▼             ▼
Data Preproc.  CellChat       CellPhoneDB
    │             Analysis       Analysis
    │             │             │
    ▼             ▼             ▼
Normalized     Interaction    Interaction
    │          Results        Results
    │             │             │
    └─────────────┴─────────────┘
                  │
                  ▼
            Build Gold Standard
                  │
                  ▼
          Generate Embeddings
                  │
                  ▼
     GeneCompass Train & Predict
                  │
                  ▼
        Cell Interaction Network
```

## Environment Setup

### 1. Python Environment

```bash
# Create conda environment
conda create -n ccc python=3.10
conda activate ccc

# Install dependencies
pip install -r requirements.txt
```

### 2. R Environment (for CellChat and CellPhoneDB)

```bash
# Install necessary R packages
R

# In R, execute:
install.packages(c("Seurat", "dplyr", "ggplot2", "patchwork"))
if (!requireNamespace("BiocManager", quietly = TRUE))
    install.packages("BiocManager")
BiocManager::install("ComplexHeatmap")

# Install CellChat
install.packages("devtools")
devtools::install_github("sqjin/CellChat")

# Install CellPhoneDB
# Download CellPhoneDB and install to conda environment
# See CellPhoneDB official documentation
```

## Data Preparation

### 1. Download Pretrained Models

Download pretrained models from the official GeneCompass repository:

- GeneCompass_Base: [Download Link](https://www.scidb.cn/en/anonymous/SUZOdk1y)
- GeneCompass_Small: [Download Link](https://www.scidb.cn/en/anonymous/SUZOdk1y)

Place model files in the `./pretrained_models/` directory.

### 2. Prepare Prior Knowledge

Download the following files from the GeneCompass repository and place them in the `./prior_knowledge/` directory:

- `human_mouse_tokens.pickle`
- `human_gene_median_after_filter.pickle`
- Other necessary knowledge files

### 3. Prepare Single-Cell Data

Ensure your single-cell data is in `.h5ad` format and contains:
- Expression matrix (`adata.X`)
- Cell metadata (`adata.obs`) with cell type column
- Gene information (`adata.var`)

## Complete Analysis Pipeline

### Step 1: Data Preprocessing (Optional)

If your data is already normalized, you can skip this step.

```bash
# Data filtering
python preprocess/filter.py

# Data normalization
python preprocess/normalized.py
```

### Step 2: Run CellChat Analysis

```bash
cd CellChatAnalysis

# 2.1 Convert h5ad to CSV
python h5ad_to_csv_fixed.py \
  --input ../your_data.h5ad \
  --output ./cellchat_output \
  --celltype_col cell_type

# 2.2 Create Seurat object
Rscript csv_to_rds_fixed.R ./cellchat_output 3 200

# 2.3 Run CellChat analysis
Rscript CellChatAnalysis_fixed.R \
  ./cellchat_output/seurat_obj.rds \
  ./cellchat_output \
  cell_type \
  20 \
  4 \
  10 \
  FALSE
```

**Output Files**:
- `cellchat_communication.csv` - Cell communication results
- `cellchat_pathways.csv` - Signaling pathway results
- `cell_interaction_strength_matrix.csv` - Interaction strength matrix

### Step 3: Run CellPhoneDB Analysis

```bash
cd CellPhoneAnalysis

# 3.1 Generate DEGs file
python prepare_DEGs_h5ad.py \
  --mode degs_only \
  --h5ad ../your_data.h5ad \
  --outdir ./cpdb_output \
  --groupby cell_type \
  --n_top_genes 250

# 3.2 Generate microenvironment file
python prepare_microenvs_h5ad.py \
  --h5ad ../your_data.h5ad \
  --output ./cpdb_output/microenv.tsv \
  --groupby cell_type

# 3.3 Run CellPhoneDB analysis
python CellPhoneAnalysis.py \
  --h5ad ../your_data.h5ad \
  --cpdb ./v5.0.0/cellphonedb.zip \
  --degs ./cpdb_output/DEGs.tsv \
  --microenv ./cpdb_output/microenv.tsv \
  --groupby cell_type \
  --outdir ./cpdb_output
```

**Output Files**:
- `significant_means.txt` - Significant interaction results
- `means.txt` - All interaction means
- `pvalues.txt` - Significance test results

### Step 4: Build Gold Standard

```bash
cd ..

python building_gold_standard_database.py \
  --cellchat_dir ./CellChatAnalysis/cellchat_output/ \
  --cpdb_dir ./CellPhoneAnalysis/cpdb_output/ \
  --output_dir ./gold_standard/ \
  --threshold_method quantile \
  --threshold_value 0.7
```

**Output Files**:
- `complete_labeled_interactions.csv` - Complete labeled matrix
- `machine_learning_dataset.csv` - Machine learning dataset
- `gold_standard_interactions.csv` - Gold standard positive samples
- `dataset_statistics.csv` - Dataset statistics
- Multiple visualization charts (`.png`)

### Step 5: Generate Embeddings

```bash
python generate_embeddings.py \
  --dataset_path ./normalized_data/TabulaSapiens/tabula_sapiens_liver/ \
  --model_path ./pretrained_models/GeneCompass_Base \
  --token_dict_path ./prior_knowledge/human_mouse_tokens.pickle \
  --output_path ./embeddings/embeddings.pickle \
  --batch_size 128 \
  --gpu_ids "0"
```

**Note**: This step may take a long time depending on the dataset size.

### Step 6: Cell-Cell Interaction Analysis

```bash
python cell_cell_interaction.py
```

**Output Files**:
- `best_model.pt` - Best trained model
- `evaluation_results.pkl` - Evaluation results
- `predictions.csv` - Prediction results (if inference code is added)

## Output File Description

### CellChat Outputs
- `cellchat_communication.csv`: Each row is a cell communication event with sender, receiver, interaction name, and probability
- `cellchat_pathways.csv`: Signaling pathways and their communication probabilities
- `cell_interaction_strength_matrix.csv`: Matrix of interaction strength between cell types
- `communication_network.png`: Cell communication network visualization
- `communication_heatmap.png`: Interaction strength heatmap

### CellPhoneDB Outputs
- `significant_means.txt`: Significant ligand-receptor interactions and their strength
- `means.txt`: Mean strength of all ligand-receptor interactions
- `pvalues.txt`: P-values from significance tests

### Gold Standard Outputs
- `complete_labeled_interactions.csv`: Complete dataset with all cell pairs and their labels
- `machine_learning_dataset.csv`: Dataset for machine learning with features and labels
- `gold_standard_interactions.csv`: Only high-confidence gold standard positive samples
- `dataset_statistics.csv`: Statistics of the dataset
- `complete_interaction_matrix.png`: Complete interaction matrix heatmap
- `gold_standard_labels.png`: Gold standard label heatmap
- `score_distributions.png`: Consensus score distribution plots
- `gold_standard_network.png`: Gold standard interaction network graph

### GeneCompass Outputs
- `embeddings.pickle`: Cell/gene embedding vectors
- `best_model.pt`: Best trained model
- `evaluation_results.pkl`: Evaluation metrics and results

## Common Issues

### Q1: Out of memory?

**A**:
- Reduce `--batch_size` parameter in `generate_embeddings.py`
- Reduce number of parallel workers in CellChat
- Use data subsampling in CellPhoneDB analysis

### Q2: Gene name format mismatch?

**A**:
- Check if your data uses Ensembl ID or Gene Symbol
- CellChat and CellPhoneDB typically use Gene Symbol
- If using Ensembl ID, convert to Gene Symbol first

### Q3: No cell-cell interactions detected?

**A**:
- Check number of cell types (at least 2 required)
- Check if each cell type has enough cells
- Lower `min_cells` parameter
- Check if gene name format matches database

### Q4: Training takes too long?

**A**:
- Use GPU acceleration
- Reduce dataset size
- Adjust model complexity
- Use GeneCompass_Small model instead of GeneCompass_Base

### Q5: How to evaluate model performance?

**A**:
- Check training and validation accuracy
- Use evaluation metrics: accuracy, precision, recall, F1-score
- Check confusion matrix
- Compare results with traditional methods (CellChat, CellPhoneDB)

## Next Steps

- Check detailed `README.md` for more features
- Check `README.md` in each subdirectory for detailed usage
- Adjust parameters according to your specific needs

## Citation

If you use this tool in your research, please cite:

- GeneCompass: [Paper Link]
- CellChat: Jin, S., et al. (2021). Nature Protocols.
- CellPhoneDB: Vento-Tormo, R., et al. (2018). Science.

## Contact

For questions or suggestions, please submit an Issue or Pull Request.
