# Technical Notes

## Core Technical Features

### GeneCompass-Based Regression Fine-Tuning

This project adopts a regression-based fine-tuning approach on the GeneCompass foundation model for cell-cell interaction analysis, rather than using simple classification models.

#### Key Technical Advantages

1. **Direct Fine-Tuning of Pre-trained Model**
   - Fine-tuned on the GeneCompass pre-trained model
   - Preserves deep feature representations learned from large-scale single-cell data
   - Avoids loss of feature expressiveness when training small models from scratch

2. **Regression Task Design**
   - Predicts continuous cell interaction strength scores
   - Uses consensus scores from gold standard as training labels
   - Better captures complexity and strength variations in cell interactions

3. **Multi-GPU Distributed Training**
   - Supports single-machine multi-GPU training
   - Uses PyTorch distributed training framework
   - Automatic data synchronization and model parallelization

#### Model Architecture

```
GeneCompass Pre-trained Model
├── BERT Encoder (Multi-layer Transformer)
│   ├── Self-attention mechanism
│   ├── Feed-forward neural network
│   └── Residual connections and layer normalization
├── Prior Knowledge Injection
│   ├── Promoter region similarity
│   ├── Co-expression correlation
│   ├── Gene family information
│   ├── Regulatory network
│   └── Homologous gene mapping
└── Regression Prediction Head
    └── Single output node (interaction strength score)
```

#### Training Pipeline

1. **Data Preparation**
   - Extract cell pairs and consensus scores from gold standard
   - Build cell pair sequences (sender sequence + receiver sequence)
   - Load pre-trained GeneCompass model and prior knowledge

2. **Model Fine-Tuning**
   - Use HuggingFace Transformers framework
   - Regression loss function: MSELoss
   - Learning rate scheduler: Linear decay with warmup
   - Evaluation metrics: MSE, RMSE, MAE, R², Correlation

3. **Prediction and Visualization**
   - Predict cell interaction scores using fine-tuned model
   - Generate interaction strength matrix
   - Draw heatmaps and distribution plots

#### Comparison with Traditional Methods

| Feature | Traditional Classification | GeneCompass Regression Fine-Tuning |
|---------|--------------------------|-----------------------------------|
| Model Basis | Simple neural network | Pre-trained large model |
| Task Type | Binary classification | Regression |
| Output | Interaction probability | Interaction strength score |
| Feature Learning | From scratch | Transfer learning |
| Prior Knowledge | None | Multiple biological priors integrated |
| Prediction Granularity | Interaction/no interaction | Continuous strength value |

#### Performance Advantages

- **Better Feature Representation**: Inherits gene representations learned by GeneCompass on large-scale data
- **More Accurate Prediction**: Regression tasks better model the continuity of interaction strength
- **Stronger Generalization**: Pre-trained model provides good initialization
- **Biological Interpretability**: Integrates multiple prior knowledge for enhanced interpretability

## Applicable Scenarios

- Cell-cell interaction analysis of single-cell RNA sequencing data
- Research requiring precise prediction of interaction strength
- Cross-species cell interaction comparison
- Efficient analysis of large-scale single-cell datasets

## Technical Requirements

- Python 3.8+
- PyTorch 1.12+
- Transformers 4.20+
- GPU (recommended, with CUDA support)
- Sufficient memory (16GB+ recommended)

## References

- GeneCompass: A foundation model for single-cell RNA sequencing analysis
- CellChat: Inferring and analyzing cell-cell communication
- CellPhoneDB: Inferring cell-cell communication from combined expression of multi-subunit ligand-receptor complexes
