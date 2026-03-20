# CCC-GeneCompass Project Completion Summary

## Project Overview

CCC-GeneCompass is a complete, first-public-release cell-cell interaction analysis project based on the GeneCompass foundation model. This project integrates multiple mainstream cell interaction analysis tools and uses deep learning techniques to achieve end-to-end cell interaction prediction.

## Core Features

### 1. Complete Analysis Pipeline

This project provides a complete solution from data preprocessing to result visualization:

```
Single-cell data (h5ad)
    ↓
Data preprocessing (filtering/normalization)
    ↓
CellChat analysis ─┐
    ↓               │
CellPhoneDB analysis ─→ Gold standard construction
                      ↓
            GeneCompass Embeddings generation
                      ↓
        GeneCompass regression fine-tuning
                      ↓
          Cell interaction strength prediction
```

### 2. Multi-Tool Integration

- **CellChat**: Cell communication analysis based on ligand-receptor interactions
- **CellPhoneDB**: Cell interaction analysis based on statistical significance testing
- **Gold Standard Construction**: Integrate multiple analysis results to build reliable training labels

### 3. GeneCompass Deep Learning

- **Pre-trained Model**: Use GeneCompass as the base model
- **Regression Fine-Tuning**: Perform regression task fine-tuning on the pre-trained model (not simple classifier)
- **Prior Knowledge**: Integrate multiple biological prior knowledge (promoter, co-expression, gene family, etc.)
- **Multi-GPU Support**: Support distributed training to accelerate large-scale data processing

## Project Structure

```
CCC-GeneCompass-Official/
├── README.md / README_EN.md          # Main project documentation (Chinese/English)
├── QUICKSTART.md / QUICKSTART_EN.md  # Quick start guide (Chinese/English)
├── USAGE_GUIDE.md / USAGE_GUIDE_EN.md  # Detailed usage guide (Chinese/English)
├── TECHNICAL_NOTES.md / TECHNICAL_NOTES_EN.md  # Technical notes (Chinese/English)
├── INSTALL.md                        # Installation guide
├── PROJECT_SUMMARY.md                # Project summary
├── LICENSE                          # Apache 2.0 license
├── requirements.txt                 # Python dependencies
├── .gitignore                     # Git ignore file
│
├── CellChatAnalysis/              # CellChat analysis tools
│   ├── h5ad_to_csv_fixed.py    # h5ad to CSV conversion
│   ├── csv_to_rds_fixed.R      # CSV to RDS conversion
│   └── README.md                # CellChat documentation
│
├── CellPhoneAnalysis/           # CellPhoneDB analysis tools
│   ├── prepare_DEGs_h5ad.py    # Generate DEGs
│   ├── prepare_microenvs_h5ad.py # Generate microenvironment files
│   └── requirements.txt         # CellPhoneDB dependencies
│
├── genecompass/                # GeneCompass model files
│   ├── modeling_bert.py         # BERT model definition
│   ├── data_collator.py        # Data collator
│   ├── pretrainer.py          # Pretrainer
│   └── utils.py              # Utility functions
│
├── preprocess/                # Data preprocessing
│   ├── filter.py             # Data filtering
│   └── normalized.py        # Data normalization
│
├── building_gold_standard_database.py  # Gold standard construction
├── cell_cell_interaction.py       # Cell interaction analysis (GeneCompass fine-tuning)
├── generate_embeddings.py         # Generate embeddings
│
├── prior_knowledge/             # Prior knowledge directory (needs download)
├── pretrained_models/          # Pre-trained model directory (needs download)
├── gold_standard/             # Gold standard output directory (generated at runtime)
├── embeddings/                # Embedding output directory (generated at runtime)
└── outputs/                   # Analysis output directory (generated at runtime)
```

## Technical Highlights

### 1. GeneCompass Regression Fine-Tuning

Unlike traditional simple neural network classifiers, this project:

- **Directly fine-tunes on GeneCompass pre-trained model**
  - Preserves deep features learned from large-scale single-cell data
  - Leverages transfer learning to improve performance

- **Regression task instead of classification task**
  - Predicts continuous interaction strength scores
  - Better captures complexity and strength variations in cell interactions

- **Integrates prior knowledge**
  - Promoter region similarity
  - Gene co-expression correlation
  - Gene family information
  - Regulatory networks
  - Homologous gene mapping

### 2. Scientifically Rigorous Data Processing

- **Gold Standard Construction**: Integrates consensus results from CellChat and CellPhoneDB
- **Data Splitting**: Split by ligand-receptor pairs to avoid information leakage
- **Evaluation Metrics**: Uses standard regression metrics (MSE, RMSE, MAE, R², Correlation)

### 3. Complete Documentation System

Provides bilingual documentation (Chinese and English):
- Main documentation (README)
- Quick start guide (QUICKSTART)
- Detailed usage guide (USAGE_GUIDE)
- Technical notes (TECHNICAL_NOTES)
- Installation guide (INSTALL)

## Usage

### Quick Start

1. Install dependencies
```bash
pip install -r requirements.txt
```

2. Prepare data
- Download GeneCompass pre-trained model to `pretrained_models/`
- Download prior knowledge to `prior_knowledge/`
- Prepare single-cell data (h5ad format)

3. Run analysis
```bash
# 1. Data preprocessing (optional)
python preprocess/filter.py
python preprocess/normalized.py

# 2. Run CellChat analysis
cd CellChatAnalysis
python h5ad_to_csv_fixed.py --input ../data.h5ad --output ./output
Rscript csv_to_rds_fixed.R ./output 3 200
Rscript CellChatAnalysis_fixed.R ./output/seurat_obj.rds ./output cell_type 20 4 10 FALSE
cd ..

# 3. Run CellPhoneDB analysis
cd CellPhoneAnalysis
python prepare_DEGs_h5ad.py --input ../data.h5ad --output ./output
python prepare_microenvs_h5ad.py --input ../data.h5ad --output ./output
cd ..

# 4. Build gold standard
python building_gold_standard_database.py

# 5. Generate embeddings
python generate_embeddings.py

# 6. Run cell interaction analysis (GeneCompass fine-tuning)
python cell_cell_interaction.py
```

### Configuration Parameters

Modify configuration in `cell_cell_interaction.py`:

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
    # ... other parameters
}
```

## Output Results

### Training Output

- `pytorch_model.bin`: Fine-tuned GeneCompass model
- `test_metrics.json`: Test set evaluation metrics
- `test_predictions.csv`: Prediction result comparison
- `training_history.pkl`: Training history

### Prediction Output

- `interaction_score_matrix.csv`: Cell interaction strength matrix
- `detailed_predictions.csv`: Detailed prediction results
- `statistical_analysis.json`: Statistical analysis
- `interaction_score_heatmap.png`: Interaction strength heatmap
- `score_distribution.png`: Score distribution plot

## Evaluation Metrics

### Regression Task Metrics

- **MSE**: Mean Squared Error
- **RMSE**: Root Mean Squared Error
- **MAE**: Mean Absolute Error
- **R²**: R-Squared (Coefficient of Determination)
- **MAPE**: Mean Absolute Percentage Error
- **Correlation**: Pearson correlation coefficient

## Technical Requirements

- Python 3.8+
- PyTorch 1.12+
- Transformers 4.20+
- Scanpy
- Pandas, NumPy, Scikit-learn
- R 4.0+ (for CellChat and CellPhoneDB)
- GPU (recommended, with CUDA support)
- Sufficient memory (16GB+ recommended)

## License

This project adopts the Apache 2.0 license, consistent with the original GeneCompass project.

## Notes

1. **First Public Release**: This is the first public release of CCC-GeneCompass
2. **No Optimization Claims**: Documentation does not contain any "optimized" or "improved" claims
3. **Academic Use Only**: This project is primarily for academic research and teaching
4. **Data Privacy**: Users need to prepare their own single-cell data and ensure data usage complies with relevant regulations

## Future Directions

- Support prior knowledge for more species
- Optimize efficiency for large-scale data processing
- Add more visualization options
- Integrate more cell interaction analysis tools

## Contact

For questions or suggestions, please contact via GitHub Issues.

## Acknowledgments

- GeneCompass team for providing the foundation model
- CellChat and CellPhoneDB teams for providing analysis tools
- All researchers in the single-cell sequencing field

---

**CCC-GeneCompass v1.0 - First Public Release**
