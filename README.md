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
