<<<<<<< HEAD
# CCC-GeneCompass

Cell-Cell Interaction Analysis Based on GeneCompass Foundation Model

## Overview

CCC-GeneCompass is a comprehensive tool for analyzing cell-cell interactions in single-cell data using the GeneCompass foundation model. This project integrates multiple state-of-the-art cell communication analysis methods to build a robust gold standard and predict cell interactions using deep learning approaches.

## Key Features

- **CellChat Integration**: Analyze cell-cell communications using ligand-receptor interactions
- **CellPhoneDB Integration**: Identify significant interactions using statistical methods
- **Gold Standard Construction**: Build consensus gold standards from multiple analysis tools
- **GeneCompass Embeddings**: Generate high-dimensional embeddings using pretrained GeneCompass model
- **GeneCompass Fine-Tuning**: Fine-tune GeneCompass model using regression tasks for cell interaction prediction
  - Direct fine-tuning on GeneCompass model, not using simple classifiers
  - Predict cell interaction strength scores using regression approach
  - Preserve GeneCompass's deep feature learning capabilities

## Analysis Pipeline

```
Single-cell Data (.h5ad)
    │
    ├─────────────┬─────────────┬─────────────┐
    ▼             ▼             ▼             ▼
Data      CellChat      CellPhoneDB
Preproc.    Analysis       Analysis
    │             │             │
    ▼             ▼             ▼
Normalized    Interaction    Interaction
    │          Results        Results
    │             │             │
    └─────────────┴─────────────┘
                  │
                  ▼
            Gold Standard
            Construction
                  │
                  ▼
          Generate Embeddings
                  │
                  ▼
        GeneCompass Training
                  │
                  ▼
        Interaction Prediction
```

## Project Structure

```
CCC-GeneCompass/
├── CellChatAnalysis/              # CellChat analysis tools
│   ├── h5ad_to_csv_fixed.py       # h5ad to CSV conversion
│   ├── csv_to_rds_fixed.R         # CSV to RDS conversion
│   ├── CellChatAnalysis_fixed.R   # CellChat main script
│   └── README.md                  # CellChat documentation
├── CellPhoneAnalysis/             # CellPhoneDB analysis tools
│   ├── prepare_DEGs_h5ad.py       # Generate DEGs
│   ├── prepare_microenvs_h5ad.py  # Generate microenvironment file
│   ├── CellPhoneAnalysis.py       # CellPhoneDB main script
│   └── requirements.txt            # CellPhoneDB dependencies
├── genecompass/                   # GeneCompass model files
│   ├── modeling_bert.py           # BERT model definition
│   ├── data_collator.py          # Data collator
│   └── utils.py                  # Utility functions
├── preprocess/                    # Data preprocessing
│   ├── filter.py                  # Data filtering
│   └── normalized.py             # Data normalization
├── gold_standard/                # Gold standard output directory
├── embeddings/                   # Embeddings output directory
├── outputs/                      # Analysis output directory
├── building_gold_standard_database.py  # Gold standard construction
├── cell_cell_interaction.py       # Main analysis script
├── generate_embeddings.py         # Embedding generation
├── requirements.txt               # Python dependencies
├── README.md                     # Main documentation (Chinese)
├── README_EN.md                  # Main documentation (English)
├── QUICKSTART.md                 # Quick start guide (Chinese)
├── QUICKSTART_EN.md              # Quick start guide (English)
└── LICENSE                       # License file
```

## Quick Start

### 1. Environment Setup

```bash
# Python environment
conda create -n ccc python=3.10
conda activate ccc
pip install -r requirements.txt

# R environment (for CellChat and CellPhoneDB)
# Install necessary R packages
```

### 2. Data Preparation

- Download pretrained GeneCompass models
- Prepare prior knowledge files
- Prepare single-cell data in `.h5ad` format

### 3. Run Analysis Pipeline

For detailed instructions, see [QUICKSTART_EN.md](QUICKSTART_EN.md).

## Core Components

### CellChat Analysis

CellChat is a powerful tool for inferring, analyzing, and visualizing cell-cell communication.

**Key Features**:
- Database of known ligand-receptor interactions
- Probabilistic model for communication inference
- Visualization tools for interaction networks

**Usage**:
```bash
cd CellChatAnalysis
python h5ad_to_csv_fixed.py --input data.h5ad --output ./output
Rscript csv_to_rds_fixed.R ./output
Rscript CellChatAnalysis_fixed.R ./output/seurat_obj.rds ./output
```

### CellPhoneDB Analysis

CellPhoneDB uses statistical methods to identify significant cell-cell interactions.

**Key Features**:
- Statistical significance testing
- Support for custom databases
- Comprehensive interaction database

**Usage**:
```bash
cd CellPhoneAnalysis
python prepare_DEGs_h5ad.py --h5ad data.h5ad --outdir ./output
python prepare_microenvs_h5ad.py --h5ad data.h5ad --output ./output/microenv.tsv
python CellPhoneAnalysis.py --h5ad data.h5ad --outdir ./output
```

### Gold Standard Construction

Integrate results from CellChat and CellPhoneDB to build a consensus gold standard.

**Method**:
- Normalize interaction scores from both tools
- Calculate consensus score: (Norm_CellChat + Norm_CPDB_Mean + Norm_CPDB_Max) / 3
- Apply threshold to identify high-confidence interactions

**Usage**:
```bash
python building_gold_standard_database.py \
  --cellchat_dir ./cellchat_output \
  --cpdb_dir ./cpdb_output \
  --output_dir ./gold_standard \
  --threshold_method quantile \
  --threshold_value 0.7
```

### GeneCompass Embeddings

Generate high-dimensional embeddings for cells using the pretrained GeneCompass model.

**Usage**:
```bash
python generate_embeddings.py \
  --dataset_path ./normalized_data/ \
  --model_path ./pretrained_models/GeneCompass_Base \
  --token_dict_path ./prior_knowledge/human_mouse_tokens.pickle \
  --output_path ./embeddings/embeddings.pickle
```

### Cell-Cell Interaction Prediction

Train a neural network model to predict cell-cell interactions using embeddings and gold standard labels.

**Usage**:
```bash
python cell_cell_interaction.py
```

## Output Files

### CellChat Outputs
- `cellchat_communication.csv`: Cell communication events
- `cellchat_pathways.csv`: Signaling pathways
- `cell_interaction_strength_matrix.csv`: Interaction strength matrix
- `communication_network.png`: Network visualization
- `communication_heatmap.png`: Heatmap visualization

### CellPhoneDB Outputs
- `significant_means.txt`: Significant interactions
- `means.txt`: All interaction means
- `pvalues.txt`: Significance test results

### Gold Standard Outputs
- `complete_labeled_interactions.csv`: Complete labeled dataset
- `machine_learning_dataset.csv`: ML dataset with features
- `gold_standard_interactions.csv`: High-confidence positive samples
- `dataset_statistics.csv`: Dataset statistics
- Various visualizations (`.png` files)

### GeneCompass Outputs
- `embeddings.pickle`: Cell/gene embedding vectors
- `best_model.pt`: Trained model
- `evaluation_results.pkl`: Evaluation metrics

### GeneCompass Pretrained-model

Pretrained models of GeneCompass on 100 million single-cell transcriptomes from humans and mice. Put pretrained_model dir under main path.('./pretrained_models/GeneCompass_Small', './pretrained_models/GeneCompass_Base')

| Model             | Description                         | Download                                           |
| ----------------- | ----------------------------------- | -------------------------------------------------- |
| GeneCompass_Small | Pretrained on 6-layer GeneCompass.  | [Link](https://www.scidb.cn/en/anonymous/SUZOdk1y) |
| GeneCompass_Base  | Pretrained on 12-layer GeneCompass. | [Link](https://www.scidb.cn/en/anonymous/SUZOdk1y) |

## Requirements

### System Requirements
- Python: 3.8+
- R: 4.0+
- Memory: At least 16GB (recommended 32GB+)
- GPU: Optional, for accelerating model training

### Python Dependencies
See `requirements.txt` for complete list.

### R Dependencies
- Seurat
- CellChat
- dplyr
- ggplot2
- ComplexHeatmap
- patchwork

## Citation

If you use this tool in your research, please cite:

- **GeneCompass**: Yang, X., Liu, G., Feng, G. *et al.* GeneCompass: deciphering universal gene regulatory mechanisms with a knowledge-informed cross-species foundation model. *Cell Res* **34**, 830–845 (2024). https://doi.org/10.1038/s41422-024-01034-y
- **CellChat**: Jin, S., et al. (2021). Inferring cell-cell communication by integrating ligand-receptor, signaling gene, and TF-target networks. Nature Protocols.
- **CellPhoneDB**: Vento-Tormo, R., et al. (2018). Single-cell reconstruction of developmental trajectories during human endometriosis. Science.

## License

This project follows the license of the original GeneCompass project.

## Acknowledgments

- GeneCompass: Foundation model for single-cell analysis
- CellChat: Cell-cell communication inference tool
- CellPhoneDB: Database of ligand-receptor interactions

## Contact

For questions, issues, or suggestions, please:
- Submit an Issue on GitHub
- Contact the development team

---

**Version**: 1.0.2
**Release Date**: March 2026
=======
# CCC-GeneCompass: Cell-Cell Communication via GeneCompass

**Large Language Model for Cell-Cell Interaction Prediction in Single-Cell Transcriptomics**

[![Python 3.11+](https://img.shields.io/badge/python-3.11+-blue.svg)](https://www.python.org/)
[![PyTorch](https://img.shields.io/badge/PyTorch-2.5.0-red.svg)](https://pytorch.org/)
[![License](https://img.shields.io/badge/License-MIT-green.svg)](LICENSE)
[中文文档](README_zh.md)

---

## Overview

CCC-GeneCompass leverages the **GeneCompass** — a BERT-based large language model pretrained on 100M+ single-cell transcriptomes — to predict cell-cell interaction (CCI) strength. The pipeline:

1. Computes **cell-type-aggregated expression profiles** from single-cell data
2. Constructs a **gold standard** from CellChat + CellPhoneDB v5 consensus
3. Fine-tunes GeneCompass to rank cell-type interaction strength
4. Evaluates via **5-fold cross-validation** with **Spearman rank correlation ρ**

## Key Features

- **Deterministic**: Cell-type mean expression eliminates random single-cell sampling noise
- **Scientifically rigorous**: Spearman ρ + bootstrap 95% CI + permutation p-value
- **Multi-GPU**: DataParallel training on up to 4 GPUs
- **End-to-end pipeline**: Preprocessing → Gold Standard → Training → CV → Visualization
- **Modular design**: Separate scripts for each step, easy to customize

## Installation

```bash
pip install -r requirements.txt

# R dependencies for CellChat
R -e 'install.packages(c("devtools","NMF","circlize","ComplexHeatmap"))'
R -e 'devtools::install_github("jinworks/CellChat")'

# CellPhoneDB v5 database
# Download cellphonedb.zip from: https://github.com/ventolab/CellphoneDB
```

## Data Preparation

### 1. Single-Cell Data (h5ad)
Requirements: `.obs` must contain a `cell_type` column; gene symbols accessible via `feature_name` or `var_names`.

### 2. Pretrained Model
Download GeneCompass_Base checkpoint and place:
```
pretrained_models/
├── pytorch_model.bin   # ~1.1GB
└── config.json
```

### 3. Knowledge Files (for tokenization)
```
prior_knowledge/
├── human_mouse_tokens.pickle
└── public/
    └── human_gene_median_after_filter.pickle
```

### 4. CellPhoneDB v5 Database
```
CellPhoneAnalysis/v5.0.0/
├── cellphonedb.zip
├── gene_input.csv
└── protein_input.csv
```

## Pipeline

```bash
# ================================================================
# Full Pipeline for One Organ
# ================================================================
# Input:  Raw single-cell h5ad (must have cell_type column)
# Output: results/ (5-fold CV metrics + models + visualizations)
# ================================================================

ORGAN=pancreas
RAW_H5AD=/path/to/original/${ORGAN}.h5ad       # original data

# ====== Step 0: Data Preprocessing ======
#  Raw h5ad → filtered.h5ad + tokenized arrow + cell-type-aggregated
python preprocess_data.py \
    --h5ad ${RAW_H5AD} \
    --output data/${ORGAN} \
    --tokens prior_knowledge/human_mouse_tokens.pickle \
    --medians prior_knowledge/public/human_gene_median_after_filter.pickle
#  Outputs: data/${ORGAN}/filtered.h5ad              (filtered data)
#           data/${ORGAN}/single_cell_dataset/       (tokenized arrow)
#           data/${ORGAN}/cell_type_aggregated/      (used by Step 4)

# ====== Step 1: CellChat Analysis ======
#  Filtered h5ad → interaction matrix + communication probabilities
python CellChatAnalysis/h5ad_to_csv.py \
    --input data/${ORGAN}/filtered.h5ad \
    --output data/${ORGAN}/cellchat/
Rscript CellChatAnalysis/csv_to_rds.R data/${ORGAN}/cellchat/ 3 200
Rscript CellChatAnalysis/cellchat_gold_standard.R data/${ORGAN}/cellchat/ 4
#  Outputs: data/${ORGAN}/cellchat/cell_interaction_strength_matrix.csv
#           data/${ORGAN}/cellchat/cellchat_communication.csv

# ====== Step 2: CellPhoneDB v5 Analysis ======
#  Raw h5ad → statistically significant interaction matrix
CPDB_ZIP=/path/to/cellphonedb.zip
CPDB_DATA=/path/to/CellPhoneAnalysis/v5.0.0/
python run_cpdb.py \
    --h5ad ${RAW_H5AD} \
    --cpdb_db ${CPDB_ZIP} \
    --cpdb_genes ${CPDB_DATA} \
    --output data/${ORGAN}/cellphonedb/
#  Output: data/${ORGAN}/cellphonedb/significant_means.txt

# ====== Step 3: Joint Gold Standard ======
#  Weights automatically learned from source significance (w ∝ mean sig LR pairs)
python genecompass_gold_standard.py \
    --cellchat data/${ORGAN}/cellchat \
    --cpdb data/${ORGAN}/cellphonedb \
    --output data/${ORGAN}/gold_standard
#  Output: data/${ORGAN}/gold_standard/complete_labeled_interactions.csv

# ====== Step 4: 5-Fold Cross-Validation ======
#  Train GeneCompass with gold standard labels, evaluate Spearman ρ
python pipeline_cv.py \
    --proj_root . \
    --gs_path data/${ORGAN}/gold_standard/complete_labeled_interactions.csv \
    --dataset data/${ORGAN}/cell_type_aggregated \
    --output results/${ORGAN}_cv \
    --organ ${ORGAN^} --epochs 30 --batch 1 --grad_accum 4
#  Outputs: results/${ORGAN}_cv/cv_summary.json + fold{1-5}/ + visualizations

# ====== Step 5: Standalone Inference (optional) ======
python pipeline_inference.py \
    --model results/${ORGAN}_cv/fold1/best_model \
    --test_set results/${ORGAN}_cv/fold1/data_splits/test \
    --token_dict prior_knowledge/human_mouse_tokens.pickle
```

## Evaluation Metrics

### Primary: Spearman Rank Correlation ρ
Measures how well the model ranks cell-type pairs by interaction strength.

| ρ range | Interpretation |
|---------|---------------|
| ρ ≥ 0.7 | EXCELLENT |
| 0.5 ≤ ρ < 0.7 | GOOD |
| 0.3 ≤ ρ < 0.5 | MODERATE |
| ρ < 0.3 | WEAK |

Each fold reports: ρ with bootstrap 95% CI + permutation p-value.
Final: mean ± std across 5 folds.

### Secondary
- **Pearson r**: Linear correlation with gold standard
- **R²**: Explained variance
- **RMSE**: Root mean squared error

## Output Structure

```
data/{organ}/
├── filtered.h5ad                       # filtered data (Step 0)
├── single_cell_dataset/                # tokenized arrow data (Step 0)
├── cell_type_aggregated/               # cell-type-aggregated data (Step 0)
├── cellchat/                           # CellChat outputs (Step 1)
│   ├── cell_interaction_strength_matrix.csv
│   └── cellchat_communication.csv
├── cellphonedb/                        # CellPhoneDB outputs (Step 2)
│   ├── significant_means.txt
│   └── statistical_analysis_*.txt
└── gold_standard/                      # joint gold standard (Step 3)
    ├── complete_labeled_interactions.csv
    └── gold_standard_stats.json

results/{organ}_cv/                     # 5-fold CV outputs (Step 4)
├── cv_summary.json                     # mean ± std across 5 folds
├── fold{1-5}/
│   ├── metrics.json                    # per-fold ρ + CI + Pearson + R²
│   ├── best_model/                     # trained model (pytorch_model.bin)
│   ├── test_true.npy                   # true labels
│   └── test_pred.npy                   # predictions
├── interaction_heatmap.png             # 300dpi visualizations
├── interaction_network.png
├── interaction_circular.png
├── interaction_bubble.png
├── interaction_flow.png
├── autocrine_scores.png
└── true_vs_predicted.png
```

## Citation

```bibtex
@software{ccc_genecompass,
  title = {CCC-GeneCompass: Cell-Cell Communication via Large Language Model},
  year = {2025},
  note = {Based on GeneCompass: A Large-Scale Pretrained Model for Single-Cell Gene Expression}
}
```

## License

MIT
>>>>>>> 0a00b04 (Initial commit: CCC-GeneCompass v3 - gene compass BERT model for cell-cell communication analysis)
