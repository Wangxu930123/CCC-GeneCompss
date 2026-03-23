# CCC-GeneCompass Usage Guide

## First-Time Setup Checklist

### 1. Environment Check

```bash
# Python Environment
python --version  # Should be 3.8+
conda info --envs

# Check dependencies
python -c "import torch; import transformers; import scanpy; print('✓ Python OK')"

# R Environment
R --version  # Should be 4.0+

# In R, check:
library(Seurat)
library(CellChat)
print('✓ R OK')
```

### 2. Data Preparation

**Required Files:**
- [ ] Pretrained models in `./pretrained_models/` directory
- [ ] Token dictionary in `./prior_knowledge/` directory
- [ ] Single-cell data (`.h5ad` format) ready
- [ ] Data contains cell type information (`cell_type` column)

### 3. Directory Structure Verification

```bash
# Run this to check directory structure
ls -la
ls CellChatAnalysis/
ls CellPhoneAnalysis/
ls genecompass/
ls preprocess/
```

## Complete Analysis Workflow

### Phase 1: Data Preprocessing

**Optional Step** - Skip if your data is already normalized

```bash
# Step 1.1: Data Filtering
cd preprocess
# Modify paths in filter.py to point to your data
python filter.py

# Step 1.2: Data Normalization
# Modify paths in normalized.py to point to your data
python normalized.py
cd ..
```

**Expected Output:**
- `./filtered_data/` - Filtered data
- `./normalized_data/` - Normalized data (contains dataset folder)

### Phase 2: CellChat Analysis

```bash
cd CellChatAnalysis

# Step 2.1: Convert h5ad to CSV
python h5ad_to_csv_fixed.py \
  --input ../path/to/your_data.h5ad \
  --output ./cellchat_output \
  --celltype_col cell_type

# Step 2.2: Convert CSV to Seurat object
Rscript csv_to_rds_fixed.R ./cellchat_output 3 200

# Step 2.3: Run CellChat analysis
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

**Expected Output:**
- `./CellChatAnalysis/cellchat_output/cellchat_communication.csv`
- `./CellChatAnalysis/cellchat_output/cellchat_pathways.csv`
- `./CellChatAnalysis/cellchat_output/cell_interaction_strength_matrix.csv`
- Multiple visualization plots

**Common Issues:**
- If R errors with "multicore not supported", script will automatically detect and use "multisession"
- If gene names are Ensembl IDs, convert to Gene Symbols first
- If no interactions detected, check if you have enough cell types (at least 2)

### Phase 3: CellPhoneDB Analysis

```bash
cd CellPhoneAnalysis

# Step 3.1: Generate DEGs
python prepare_DEGs_h5ad.py \
  --mode degs_only \
  --h5ad ../path/to/your_data.h5ad \
  --outdir ./cpdb_output \
  --groupby cell_type \
  --n_top_genes 250

# Step 3.2: Generate microenv file
python prepare_microenvs_h5ad.py \
  --h5ad ../path/to/your_data.h5ad \
  --output ./cpdb_output/microenv.tsv \
  --groupby cell_type

# Step 3.3: Run CellPhoneDB
python CellPhoneAnalysis.py \
  --h5ad ../path/to/your_data.h5ad \
  --cpdb ./path/to/cellphonedb.zip \
  --degs ./cpdb_output/DEGs.tsv \
  --microenv ./cpdb_output/microenv.tsv \
  --groupby cell_type \
  --outdir ./cpdb_output

cd ..
```

**Expected Output:**
- `./CellPhoneAnalysis/cpdb_output/significant_means.txt`
- `./CellPhoneAnalysis/cpdb_output/means.txt`
- `./CellPhoneAnalysis/cpdb_output/pvalues.txt`

**Note:**
- Need to install CellPhoneDB first, refer to CellPhoneDB official documentation
- Ensure CellPhoneDB database file path is correct

### Phase 4: Build Gold Standard

```bash
# Step 4: Build gold standard dataset
python building_gold_standard_database.py \
  --cellchat_dir ./CellChatAnalysis/cellchat_output/ \
  --cpdb_dir ./CellPhoneAnalysis/cpdb_output/ \
  --output_dir ./gold_standard/ \
  --threshold_method quantile \
  --threshold_value 0.7
```

**Expected Output:**
- `./gold_standard/complete_labeled_interactions.csv`
- `./gold_standard/machine_learning_dataset.csv`
- `./gold_standard/gold_standard_interactions.csv`
- `./gold_standard/dataset_statistics.csv`
- Multiple visualization charts

**Parameter Explanation:**
- `--threshold_method quantile`: Use quantile threshold (recommended)
- `--threshold_value 0.7`: Top 30% of interactions as gold standard positives
- Can adjust to other values (e.g., 0.5 or 0.8) to change number of positive samples

### Phase 5: Generate Embeddings

```bash
# Step 5: Generate GeneCompass embeddings
python generate_embeddings.py \
  --dataset_path ./normalized_data/TabulaSapiens/tabula_sapiens_liver/ \
  --model_path ./pretrained_models/GeneCompass_Base \
  --token_dict_path ./prior_knowledge/human_mouse_tokens.pickle \
  --output_path ./embeddings/embeddings.pickle \
  --batch_size 128 \
  --gpu_ids "0"
```

**Expected Output:**
- `./embeddings/embeddings.pickle` - Embedding vectors file

**Note:**
- This step may take a long time (several hours), depending on dataset size
- If out of memory, reduce `--batch_size`
- If no GPU, omit `--gpu_ids` parameter or set to empty string

### Phase 6: Cell-Cell Interaction Analysis (GeneCompass Fine-Tuning)

```bash
# Step 6: Fine-tune GeneCompass model for regression task
# Need to modify paths in cell_cell_interaction.py
python cell_cell_interaction.py
```

**Expected Output:**
- `./outputs/pytorch_model.bin` - Fine-tuned GeneCompass model
- `./outputs/test_metrics.json` - Test evaluation metrics
- `./outputs/test_predictions.csv` - Prediction results
- `./outputs/predictions/` - Final interaction predictions and visualizations
- Training logs showing performance for each epoch

**Configuration:**
Modify following variables in `cell_cell_interaction.py`:
```python
embeddings_path = './embeddings/embeddings.pickle'
gold_standard_path = './gold_standard/machine_learning_dataset.csv'
dataset_path = './normalized_data/TabulaSapiens/tabula_sapiens_liver/'
token_dict_path = './prior_knowledge/human_mouse_tokens.pickle'
model_path = './pretrained_models/GeneCompass_Base'
output_dir = './outputs/'
```

**Important Notes:**
- This step performs regression task fine-tuning on the GeneCompass model
- Unlike simple classifiers, this directly fine-tunes the pre-trained GeneCompass model
- Predicts continuous interaction strength scores rather than binary labels
- Preserves GeneCompass's deep feature learning capabilities

## Result Interpretation

### Gold Standard Dataset

**File Structure:**
```csv
Sender,Receiver,Pair_ID,Gold_Standard_Label,Confidence_Level,Consensus_Score,...
CellA,CellB,CellA_CellB,1,Gold Standard,0.85,...
```

**Column Description:**
- `Sender`: Sender cell type
- `Receiver`: Receiver cell type
- `Gold_Standard_Label`: 1=gold standard positive, 0=negative
- `Confidence_Level`: Confidence level
- `Consensus_Score`: Consensus score (0-1)

### Evaluation Metrics (Regression Task)

- **MSE (Mean Squared Error)**: Average of squared errors
- **RMSE (Root Mean Squared Error)**: Square root of MSE
- **MAE (Mean Absolute Error)**: Average of absolute errors
- **R² (R-Squared)**: Coefficient of determination
- **MAPE (Mean Absolute Percentage Error)**: Average absolute percentage error
- **Correlation**: Pearson correlation coefficient

**Key Points:**
- This project uses GeneCompass-based regression fine-tuning, not simple classification
- The model is fine-tuned directly on the GeneCompass pre-trained model
- Prediction output is continuous interaction strength scores, not binary labels

- **Accuracy**: Proportion of correct predictions
- **Precision**: Proportion of positive predictions that are actually positive
- **Recall**: Proportion of actual positives that are predicted positive
- **F1-Score**: Harmonic mean of precision and recall

## Common Issues Troubleshooting

### Issue 1: Out of Memory

**Symptoms:**
- Python: `MemoryError` or program crashes
- R: Out of memory error

**Solutions:**
```bash
# Python: Reduce batch size
python generate_embeddings.py --batch_size 64  # Reduce from 128 to 64

# R: Reduce parallel workers
Rscript CellChatAnalysis_fixed.R ... 2  # Reduce from 4 to 2

# Close other applications
```

### Issue 2: GPU Not Used

**Symptoms:**
- Training is very slow
- `nvidia-smi` shows 0% GPU utilization

**Solutions:**
```bash
# Check CUDA availability
python -c "import torch; print(torch.cuda.is_available())"

# Check PyTorch version
python -c "import torch; print(torch.__version__)"

# If returns False, need CUDA-enabled PyTorch version
```

### Issue 3: Gene Name Mismatch

**Symptoms:**
- CellChat or CellPhoneDB detected no interactions
- "0 features" or similar error

**Solutions:**
```bash
# Check gene name format
python -c "import scanpy as sc; adata=sc.read_h5ad('data.h5ad'); print(adata.var_names[:10])"

# If Ensembl ID (like ENSG00000000003), need to convert
# Use tools in CellChatAnalysis directory
```

### Issue 4: Path Errors

**Symptoms:**
- `FileNotFoundError` or file does not exist error

**Solutions:**
```bash
# Use absolute paths
python script.py --input "/full/path/to/file.h5ad"

# Or set working directory in code
import os
os.chdir('/path/to/project')
```

## Advanced Usage

### Custom Thresholds

```bash
python building_gold_standard_database.py \
  --threshold_method absolute \
  --threshold_value 0.5  # Use absolute threshold
```

### Use Different Models

```bash
python generate_embeddings.py \
  --model_path ./pretrained_models/GeneCompass_Small  # Use small model
```

### Adjust Model Architecture

Modify model definition in `cell_cell_interaction.py`:

```python
class SimpleCellInteractionClassifier(nn.Module):
    def __init__(self, input_dim, hidden_dim=512, num_classes=2):  # Increase hidden dimension
        ...
```

## Next Steps

- Try running the complete pipeline on your own dataset
- Adjust parameters (thresholds, model architecture, etc.) based on results
- Explore different visualization options
- Compare results with traditional methods

## Getting Help

1. Check [README_EN.md](README_EN.md) for project overview
2. Check [QUICKSTART_EN.md](QUICKSTART_EN.md) for quick start
3. Check [INSTALL.md](INSTALL.md) for installation guide
4. Submit an Issue for technical support

---

**Version**: 1.0
**Last Updated**: March 2026
