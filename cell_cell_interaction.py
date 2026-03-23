#!/usr/bin/env python3
"""
Cell-Cell Interaction Analysis Based on GeneCompass

This script performs cell-cell interaction analysis using GeneCompass fine-tuning
for regression tasks, with gold standard labels derived from CellChat and CellPhoneDB.
"""

import os
import pickle
import pandas as pd
import numpy as np
from datasets import Dataset, load_from_disk
from transformers import Trainer, TrainingArguments
from genecompass import BertForSequenceClassification, DataCollatorForCellClassification
from genecompass.utils import load_prior_embedding
from sklearn.metrics import mean_squared_error, mean_absolute_error, r2_score
from sklearn.model_selection import train_test_split
import torch
import torch.distributed as dist
import logging
from collections import Counter
from tqdm import tqdm
import warnings
import matplotlib.pyplot as plt
import seaborn as sns
import json
from sklearn.metrics import confusion_matrix, classification_report
import scipy.stats as stats
from sklearn.utils.class_weight import compute_class_weight

warnings.filterwarnings('ignore')

# Configure logging
logging.basicConfig(level=logging.INFO, format='%(asctime)s - %(levelname)s - %(message)s')
logger = logging.getLogger(__name__)


def compute_regression_metrics(pred):
    """Compute regression task evaluation metrics"""
    labels = pred.label_ids
    preds = pred.predictions

    # Ensure preds is 1D array
    if len(preds.shape) > 1:
        preds = preds.flatten()

    # Ensure labels is 1D array
    if len(labels.shape) > 1:
        labels = labels.flatten()

    mse = mean_squared_error(labels, preds)
    mae = mean_absolute_error(labels, preds)
    r2 = r2_score(labels, preds)

    # Compute relative error
    abs_errors = np.abs(labels - preds)
    relative_errors = abs_errors / (np.abs(labels) + 1e-8)  # Avoid division by zero
    mape = np.mean(relative_errors) * 100  # Mean absolute percentage error

    # Compute correlation coefficient
    correlation = np.corrcoef(labels, preds)[0, 1] if len(labels) > 1 else 0

    return {
        'mse': mse,
        'rmse': np.sqrt(mse),
        'mae': mae,
        'r2': r2,
        'mape': mape,
        'correlation': correlation
    }


class RegressionCellInteractionDataset:
    """Regression task cell-cell interaction dataset construction class"""

    def __init__(self, embeddings_path, gold_standard_path, dataset_path, token_dict_path):
        """
        Initialize
        """
        self.embeddings_path = embeddings_path
        self.gold_standard_path = gold_standard_path
        self.dataset_path = dataset_path
        self.token_dict_path = token_dict_path

        # Load data
        self.load_data()

    def load_data(self):
        """Load all necessary data"""
        logger.info("Loading gene embeddings data...")
        with open(self.embeddings_path, 'rb') as f:
            self.gene_embeddings = pickle.load(f)

        logger.info("Loading gold standard labels...")
        self.gold_standard = pd.read_csv(self.gold_standard_path)

        logger.info("Loading original dataset...")
        self.original_dataset = load_from_disk(self.dataset_path)

        logger.info("Loading token dictionary...")
        with open(self.token_dict_path, 'rb') as f:
            self.token_dictionary = pickle.load(f)

        # Get cell type information
        if 'cell_type' in self.original_dataset.column_names:
            self.cell_types = self.original_dataset['cell_type']
        else:
            self.cell_types = [f"Cell_{i}" for i in range(len(self.original_dataset))]
            logger.warning("Cell type information not found, using default naming")

        # Create cell type to index mapping
        self.cell_to_indices = {}
        for idx, cell_type in enumerate(self.cell_types):
            if cell_type not in self.cell_to_indices:
                self.cell_to_indices[cell_type] = []
            self.cell_to_indices[cell_type].append(int(idx))

    def ensure_int_length(self, length_value):
        """Ensure length value is integer type"""
        if isinstance(length_value, (list, np.ndarray)):
            if len(length_value) > 0:
                return int(length_value[0])
            else:
                return 0
        elif isinstance(length_value, (int, np.integer)):
            return int(length_value)
        else:
            try:
                return int(length_value)
            except (ValueError, TypeError):
                logger.warning(f"Cannot convert length value {length_value} to integer, using default value 0")
                return 0

    def _build_cell_pair_sequence(self, dataset, sender_idx, receiver_idx, max_sequence_length):
        """Build cell pair sequence - fix empty sequence problem"""
        try:
            # Get cell data
            sender_data = dataset[sender_idx]
            receiver_data = dataset[receiver_idx]

            # Get special tokens
            cls_token = self.token_dictionary.get("<cls>", 1)
            sep_token = self.token_dictionary.get("<sep>", 2)
            pad_token = self.token_dictionary.get("<pad>", 0)

            # Safely get sequence data
            sender_input_ids = sender_data.get('input_ids', [])
            sender_values = sender_data.get('values', [])
            sender_length_raw = sender_data.get('length', 0)

            receiver_input_ids = receiver_data.get('input_ids', [])
            receiver_values = receiver_data.get('values', [])
            receiver_length_raw = receiver_data.get('length', 0)

            # Ensure sequence data is valid
            if (not sender_input_ids or not receiver_input_ids or
                    not sender_values or not receiver_values):
                return None

            # Ensure length values are integers
            sender_length = self.ensure_int_length(sender_length_raw)
            receiver_length = self.ensure_int_length(receiver_length_raw)

            # Check length validity
            if sender_length <= 0 or receiver_length <= 0:
                return None

            # Calculate available length
            available_length = max_sequence_length - 3
            total_length = sender_length + receiver_length

            if total_length > available_length:
                # Allocate length proportionally
                sender_ratio = sender_length / total_length
                sender_alloc = max(1, int(available_length * sender_ratio))
                receiver_alloc = max(1, available_length - sender_alloc)
            else:
                sender_alloc = sender_length
                receiver_alloc = receiver_length

            # Ensure allocation length is valid
            sender_alloc = max(1, min(sender_alloc, len(sender_input_ids)))
            receiver_alloc = max(1, min(receiver_alloc, len(receiver_input_ids)))

            # Truncate sequences
            sender_input_ids_trunc = sender_input_ids[:sender_alloc]
            sender_values_trunc = sender_values[:sender_alloc]
            receiver_input_ids_trunc = receiver_input_ids[:receiver_alloc]
            receiver_values_trunc = receiver_values[:receiver_alloc]

            # Check if truncated sequences are empty
            if not sender_input_ids_trunc or not receiver_input_ids_trunc:
                return None

            # Build cell pair sequence
            pair_input_ids = [cls_token]
            pair_input_ids.extend(sender_input_ids_trunc)
            pair_input_ids.append(sep_token)
            pair_input_ids.extend(receiver_input_ids_trunc)
            pair_input_ids.append(sep_token)

            # Build value sequence
            pair_values = [0.0]
            pair_values.extend(sender_values_trunc)
            pair_values.append(0.0)
            pair_values.extend(receiver_values_trunc)
            pair_values.append(0.0)

            # Pad to fixed length
            current_length = len(pair_input_ids)
            if current_length < max_sequence_length:
                pad_length = max_sequence_length - current_length
                pair_input_ids.extend([pad_token] * pad_length)
                pair_values.extend([0.0] * pad_length)
            else:
                pair_input_ids = pair_input_ids[:max_sequence_length]
                pair_values = pair_values[:max_sequence_length]

            # Final verification
            if len(pair_input_ids) != max_sequence_length or len(pair_values) != max_sequence_length:
                return None

            sequence_length = len(sender_input_ids_trunc) + len(receiver_input_ids_trunc) + 3

            return {
                'input_ids': pair_input_ids,
                'values': pair_values,
                'length': sequence_length
            }

        except Exception as e:
            logger.error(f"Failed to build cell pair sequence: {str(e)}")
            return None

    def create_cell_pair_sequences(self, max_sequence_length=2048, balance_dataset=True):
        """
        Create cell pair sequence data - regression task version
        """
        logger.info("Creating cell pair sequences (regression task)...")

        sequences = []
        labels = []
        cell_pairs = []

        # Find continuous score column - regression task uses continuous scores
        score_column = None
        for col in ['Consensus_Score', 'Interaction_Score', 'Score', 'score', 'value']:
            if col in self.gold_standard.columns:
                score_column = col
                logger.info(f"Using continuous score column: {score_column}")
                break

        if score_column is None:
            logger.error("Continuous score column not found, please check gold standard file format")
            # Try to view all columns
            logger.error(f"Available columns: {list(self.gold_standard.columns)}")
            return sequences, labels, cell_pairs

        # Adapt sender and receiver column names
        sender_column = None
        receiver_column = None

        for col in ['Sender', 'sender', 'source', 'from']:
            if col in self.gold_standard.columns:
                sender_column = col
                break

        for col in ['Receiver', 'receiver', 'target', 'to']:
            if col in self.gold_standard.columns:
                receiver_column = col
                break

        if sender_column is None or receiver_column is None:
            logger.error("Sender or receiver column not found")
            logger.error(f"Available columns: {list(self.gold_standard.columns)}")
            return sequences, labels, cell_pairs

        logger.info(f"Using sender column: {sender_column}, receiver column: {receiver_column}")

        # Filter invalid scores
        valid_data = self.gold_standard.dropna(subset=[score_column])
        valid_data = valid_data[valid_data[score_column] >= 0]  # Ensure score is non-negative

        logger.info(f"Valid data count: {len(valid_data)}")
        logger.info(f"Score statistics - min: {valid_data[score_column].min():.4f}, "
                    f"max: {valid_data[score_column].max():.4f}, "
                    f"mean: {valid_data[score_column].mean():.4f}")

        valid_sequences_count = 0
        invalid_sequences_count = 0

        for _, row in tqdm(valid_data.iterrows(), total=len(valid_data), desc="Processing cell pairs"):
            sender_type = row[sender_column]
            receiver_type = row[receiver_column]
            score = float(row[score_column])

            if sender_type in self.cell_to_indices and receiver_type in self.cell_to_indices:
                sender_idx = int(np.random.choice(self.cell_to_indices[sender_type]))
                receiver_idx = int(np.random.choice(self.cell_to_indices[receiver_type]))

                # Build sequence
                sequence = self._build_cell_pair_sequence(
                    self.original_dataset, sender_idx, receiver_idx, max_sequence_length
                )
                if sequence is not None:
                    sequences.append(sequence)
                    labels.append(score)
                    cell_pairs.append(f"{sender_type}_{receiver_type}")
                    valid_sequences_count += 1
                else:
                    invalid_sequences_count += 1
            else:
                invalid_sequences_count += 1

        if invalid_sequences_count > 0:
            logger.warning(f"Skipped {invalid_sequences_count} invalid sequences")

        logger.info(f"Created {len(sequences)} valid cell pair sequences")
        logger.info(f"Score range: {min(labels):.4f} - {max(labels):.4f}, mean: {np.mean(labels):.4f}")
        return sequences, labels, cell_pairs

    def create_huggingface_dataset(self, sequences, labels, test_size=0.2, validation_size=0.1):
        """Create HuggingFace dataset"""
        logger.info("Creating HuggingFace dataset (regression task)...")

        if len(sequences) == 0:
            logger.error("No valid sequence data")
            return None, None, None

        # Add species information (default to 0)
        species = [0] * len(sequences)

        dataset_dict = {
            'input_ids': [seq['input_ids'] for seq in sequences],
            'values': [seq['values'] for seq in sequences],
            'length': [seq['length'] for seq in sequences],
            'species': species,
            'label': labels  # Regression task uses continuous scores as labels
        }

        # Create dataset
        dataset = Dataset.from_dict(dataset_dict)

        # Split dataset
        if len(dataset) > 1:
            train_val_test = dataset.train_test_split(test_size=test_size + validation_size, seed=42)

            if len(train_val_test['test']) > 1:
                # Calculate validation and test set proportions
                val_ratio = validation_size / (test_size + validation_size)
                val_test = train_val_test['test'].train_test_split(
                    test_size=val_ratio, seed=42
                )
                train_dataset = train_val_test['train']
                val_dataset = val_test['train']
                test_dataset = val_test['test']
            else:
                train_dataset = train_val_test['train']
                val_dataset = train_val_test['test']
                test_dataset = train_val_test['test']
        else:
            logger.warning("Dataset too small, using all data as training set")
            train_dataset = dataset
            val_dataset = dataset
            test_dataset = dataset

        logger.info(f"Training set size: {len(train_dataset)}")
        logger.info(f"Validation set size: {len(val_dataset)}")
        logger.info(f"Test set size: {len(test_dataset)}")
        logger.info(f"Training set score range: {min(train_dataset['label']):.4f} - {max(train_dataset['label']):.4f}")

        return train_dataset, val_dataset, test_dataset


def setup_multigpu_training():
    """Setup multi-GPU training environment"""
    if torch.cuda.device_count() > 1:
        logger.info(f"Detected {torch.cuda.device_count()} GPUs available")

        if not dist.is_initialized():
            if 'RANK' in os.environ and 'WORLD_SIZE' in os.environ:
                rank = int(os.environ['RANK'])
                world_size = int(os.environ['WORLD_SIZE'])
                dist.init_process_group(backend='nccl', init_method='env://')
                logger.info(f"Initialized distributed training: rank={rank}, world_size={world_size}")
                return True
            else:
                logger.info("Single machine multi-GPU mode, but distributed training not initialized")
                return False
        return True
    else:
        logger.info("Single GPU training mode")
        return False


def fine_tune_regression_model(config):
    """
    Fine-tune cell interaction regression model
    """
    logger.info("Starting cell-cell interaction regression model fine-tuning...")

    # Setup multi-GPU training
    is_multigpu = setup_multigpu_training()

    # Get current process rank
    if is_multigpu and dist.is_initialized():
        rank = dist.get_rank()
        world_size = dist.get_world_size()
        logger.info(f"Current process rank: {rank}, total processes: {world_size}")
    else:
        rank = 0
        world_size = 1

    # Only main process creates output directory
    if rank == 0:
        os.makedirs(config['output_dir'], exist_ok=True)
        logger.info(f"Created output directory: {config['output_dir']}")

    # 1. Create dataset
    if rank == 0:
        logger.info("Main process creating dataset...")
        try:
            dataset_builder = RegressionCellInteractionDataset(
                embeddings_path=config['embeddings_path'],
                gold_standard_path=config['gold_standard_path'],
                dataset_path=config['dataset_path'],
                token_dict_path=config['token_dict_path']
            )

            # Create cell pair sequences
            sequences, labels, cell_pairs = dataset_builder.create_cell_pair_sequences(
                max_sequence_length=config.get('max_sequence_length', 2048),
                balance_dataset=config.get('balance_dataset', False)  # Regression task usually doesn't need balancing
            )

            # Check if dataset is valid
            if len(sequences) == 0:
                logger.error("Created dataset is empty, cannot train")
                return None, None, None

            # Create HuggingFace dataset
            train_dataset, val_dataset, test_dataset = dataset_builder.create_huggingface_dataset(
                sequences, labels,
                test_size=config.get('test_size', 0.2),
                validation_size=config.get('validation_size', 0.1)
            )

            # Save dataset
            dataset_save_path = os.path.join(config['output_dir'], "temp_dataset")
            os.makedirs(dataset_save_path, exist_ok=True)

            train_dataset.save_to_disk(os.path.join(dataset_save_path, "train"))
            val_dataset.save_to_disk(os.path.join(dataset_save_path, "val"))
            test_dataset.save_to_disk(os.path.join(dataset_save_path, "test"))
            logger.info("Dataset saved")
        except Exception as e:
            logger.error(f"Failed to create dataset: {str(e)}")
            return None, None, None
    else:
        dataset_save_path = os.path.join(config['output_dir'], "temp_dataset")
        logger.info(f"Process {rank} waiting for dataset...")

    # Synchronize all processes
    if is_multigpu and dist.is_initialized():
        dist.barrier()

    # All processes load dataset
    if rank != 0 or (rank == 0 and 'train_dataset' not in locals()):
        logger.info(f"Process {rank} loading dataset...")
        try:
            train_dataset = load_from_disk(os.path.join(dataset_save_path, "train"))
            val_dataset = load_from_disk(os.path.join(dataset_save_path, "val"))
            test_dataset = load_from_disk(os.path.join(dataset_save_path, "test"))
            logger.info(f"Process {rank} dataset loading complete")
        except Exception as e:
            logger.error(f"Process {rank} dataset loading failed: {str(e)}")
            return None, None, None

    # 2. Load prior knowledge
    logger.info("Loading prior knowledge...")
    knowledges = {}
    try:
        out = load_prior_embedding(token_dictionary_or_path=config['token_dict_path'])
        knowledges['promoter'] = out[0] if len(out) > 0 else None
        knowledges['co_exp'] = out[1] if len(out) > 1 else None
        knowledges['gene_family'] = out[2] if len(out) > 2 else None
        knowledges['peca_grn'] = out[3] if len(out) > 3 else None
        knowledges['homologous_gene_human2mouse'] = out[4] if len(out) > 4 else None
        logger.info("Prior knowledge loaded successfully")
    except Exception as e:
        logger.warning(f"Failed to load prior knowledge: {str(e)}")
        knowledges = {
            'promoter': None, 'co_exp': None, 'gene_family': None,
            'peca_grn': None, 'homologous_gene_human2mouse': None
        }

    # 3. Load pre-trained model - regression task uses num_labels=1
    logger.info("Loading pre-trained model (regression task)...")
    try:
        # Regression task, set num_labels=1
        model = BertForSequenceClassification.from_pretrained(
            config['model_path'],
            num_labels=1,  # Regression task
            output_attentions=False,
            output_hidden_states=False,
            knowledges=knowledges,
        )
        logger.info("Pre-trained model loaded successfully (regression task)")
    except Exception as e:
        logger.error(f"Failed to load pre-trained model: {str(e)}")
        return None, None, None

    # 4. Freeze some layers (optional)
    if config.get('freeze_layers', 0) > 0:
        logger.info(f"Freezing first {config['freeze_layers']} layers")
        freeze_layers = config['freeze_layers']
        if hasattr(model, 'bert') and hasattr(model.bert, 'encoder'):
            modules_to_freeze = model.bert.encoder.layer[:freeze_layers]
            for module in modules_to_freeze:
                for param in module.parameters():
                    param.requires_grad = False

    # 5. Setup training parameters - regression task
    per_device_batch_size = config.get('batch_size', 4)

    # Adjust batch size based on dataset size
    if len(train_dataset) > 100:
        per_device_batch_size = min(per_device_batch_size, 8)
    else:
        per_device_batch_size = min(per_device_batch_size, 2)

    if is_multigpu and world_size > 1:
        effective_batch_size = per_device_batch_size * world_size
        logger.info(f"Multi-GPU training: per-device batch size={per_device_batch_size}, effective batch size={effective_batch_size}")
    else:
        effective_batch_size = per_device_batch_size
        logger.info(f"Single GPU training: batch size={effective_batch_size}")

    # Adjust training epochs based on dataset size
    num_epochs = config.get('num_epochs', 40)
    if len(train_dataset) > 100:
        num_epochs = min(num_epochs, 30)

    training_args = TrainingArguments(
        output_dir=config['output_dir'],
        num_train_epochs=num_epochs,
        per_device_train_batch_size=per_device_batch_size,
        per_device_eval_batch_size=per_device_batch_size,
        learning_rate=config.get('learning_rate', 5e-5),
        evaluation_strategy="epoch",
        save_strategy="epoch",
        logging_dir=os.path.join(config['output_dir'], "logs"),
        logging_steps=10,
        disable_tqdm=False,
        lr_scheduler_type="linear",
        warmup_steps=config.get('warmup_steps', 100),
        weight_decay=config.get('weight_decay', 0.001),
        load_best_model_at_end=True,
        metric_for_best_model=config.get('metric_for_best_model', 'rmse'),
        greater_is_better=False,
        dataloader_num_workers=min(2, os.cpu_count() // 2),
        dataloader_pin_memory=True,
        fp16=config.get('fp16', True),
        local_rank=rank,
        ddp_find_unused_parameters=False,
        report_to=[],
        remove_unused_columns=True  # Key fix: automatically remove unused columns to prevent length parameter errors
    )

    # 6. Create trainer - use regression evaluation metrics
    try:
        trainer = Trainer(
            model=model,
            args=training_args,
            data_collator=DataCollatorForCellClassification(),
            train_dataset=train_dataset,
            eval_dataset=val_dataset,
            compute_metrics=compute_regression_metrics
        )
    except Exception as e:
        logger.error(f"Failed to create trainer: {str(e)}")
        return None, None, None

    # 7. Start training
    logger.info("Starting regression model training...")
    try:
        train_result = trainer.train()
        logger.info("Regression model training complete")
    except Exception as e:
        logger.error(f"Error during training: {str(e)}")
        return None, None, None

    # 8. Save model and evaluation results
    if rank == 0:
        logger.info("Saving regression model and evaluation results...")

        try:
            # Evaluate on test set
            test_predictions = trainer.predict(test_dataset)
            val_predictions = trainer.predict(val_dataset)

            # Save evaluation metrics
            with open(os.path.join(config['output_dir'], "test_metrics.json"), 'w') as f:
                json.dump(test_predictions.metrics, f, indent=2)

            # Save model
            trainer.save_model(config['output_dir'])

            # Save training history
            training_history = {
                'train_loss': [log for log in trainer.state.log_history if 'loss' in log],
                'eval_metrics': [log for log in trainer.state.log_history if 'eval_loss' in log]
            }
            with open(os.path.join(config['output_dir'], "training_history.pkl"), 'wb') as f:
                pickle.dump(training_history, f)

            # Save prediction results
            predictions_df = pd.DataFrame({
                'true_labels': test_predictions.label_ids.flatten(),
                'predicted_scores': test_predictions.predictions.flatten()
            })
            predictions_df.to_csv(os.path.join(config['output_dir'], "test_predictions.csv"), index=False)

            # Clean up temporary dataset
            if os.path.exists(dataset_save_path):
                import shutil
                shutil.rmtree(dataset_save_path)
                logger.info("Cleaned up temporary dataset")

            logger.info(f"Regression model training complete! Model and results saved to: {config['output_dir']}")
        except Exception as e:
            logger.error(f"Error saving results: {str(e)}")
            test_predictions = None
    else:
        test_predictions = None

    # Synchronize all processes
    if is_multigpu and dist.is_initialized():
        dist.barrier()

    return trainer, test_predictions, val_predictions


class RegressionCellInteractionPredictor:
    """Regression task cell-cell interaction predictor"""

    def __init__(self, model_path, token_dict_path):
        """
        Initialize predictor
        """
        self.model_path = model_path
        self.token_dict_path = token_dict_path
        self.model = None
        self.token_dictionary = None
        self.device = torch.device("cuda" if torch.cuda.is_available() else "cpu")

        # Load model and token dictionary
        self.load_model_and_tokenizer()

    def load_model_and_tokenizer(self):
        """Load model and token dictionary"""
        logger.info("Loading regression model and token dictionary...")

        # Load token dictionary
        with open(self.token_dict_path, 'rb') as f:
            self.token_dictionary = pickle.load(f)

        # Load prior knowledge
        knowledges = {}
        try:
            out = load_prior_embedding(token_dictionary_or_path=self.token_dict_path)
            knowledges['promoter'] = out[0] if len(out) > 0 else None
            knowledges['co_exp'] = out[1] if len(out) > 1 else None
            knowledges['gene_family'] = out[2] if len(out) > 2 else None
            knowledges['peca_grn'] = out[3] if len(out) > 3 else None
            knowledges['homologous_gene_human2mouse'] = out[4] if len(out) > 4 else None
        except Exception as e:
            logger.warning(f"Failed to load prior knowledge: {str(e)}")
            knowledges = {
                'promoter': None, 'co_exp': None, 'gene_family': None,
                'peca_grn': None, 'homologous_gene_human2mouse': None
            }

        # Load model - regression task
        try:
            self.model = BertForSequenceClassification.from_pretrained(
                self.model_path,
                knowledges=knowledges,
            )
            logger.info("Regression model loaded successfully")
        except Exception as e:
            logger.error(f"Failed to load regression model: {str(e)}")
            raise

        self.model.eval()
        self.model = self.model.to(self.device)
        logger.info("Regression model loading complete")

    def _create_cell_index_mapping(self, dataset):
        """Create cell type to index mapping"""
        cell_to_indices = {}
        if 'cell_type' in dataset.column_names:
            cell_types = dataset['cell_type']
        else:
            cell_types = [f"Cell_{i}" for i in range(len(dataset))]
        for idx, cell_type in enumerate(cell_types):
            if cell_type not in cell_to_indices:
                cell_to_indices[cell_type] = []
            cell_to_indices[cell_type].append(idx)
        return cell_to_indices

    def _build_cell_pair_sequence(self, dataset, sender, receiver, max_sequence_length):
        """Build cell pair sequence - complete implementation"""
        try:
            # Create cell type to index mapping
            cell_to_indices = self._create_cell_index_mapping(dataset)

            # Check if cell types exist
            if sender not in cell_to_indices or receiver not in cell_to_indices:
                logger.warning(f"Cell type {sender} or {receiver} not in dataset")
                return None

            # Randomly select indices
            sender_idx = int(np.random.choice(cell_to_indices[sender]))
            receiver_idx = int(np.random.choice(cell_to_indices[receiver]))

            # Get cell data
            sender_data = dataset[sender_idx]
            receiver_data = dataset[receiver_idx]

            # Get special tokens
            cls_token = self.token_dictionary.get("<cls>", 1)
            sep_token = self.token_dictionary.get("<sep>", 2)
            pad_token = self.token_dictionary.get("<pad>", 0)

            # Safely get sequence data
            sender_input_ids = sender_data.get('input_ids', [])
            sender_values = sender_data.get('values', [])
            sender_length = len(sender_input_ids)

            receiver_input_ids = receiver_data.get('input_ids', [])
            receiver_values = receiver_data.get('values', [])
            receiver_length = len(receiver_input_ids)

            # Ensure sequence data is valid
            if (not sender_input_ids or not receiver_input_ids or
                    not sender_values or not receiver_values):
                return None

            # Calculate available length
            available_length = max_sequence_length - 3
            total_length = sender_length + receiver_length

            if total_length > available_length:
                # Allocate length proportionally
                sender_ratio = sender_length / total_length
                sender_alloc = max(1, int(available_length * sender_ratio))
                receiver_alloc = max(1, available_length - sender_alloc)
            else:
                sender_alloc = sender_length
                receiver_alloc = receiver_length

            # Ensure allocation length is valid
            sender_alloc = max(1, min(sender_alloc, len(sender_input_ids)))
            receiver_alloc = max(1, min(receiver_alloc, len(receiver_input_ids)))

            # Truncate sequences
            sender_input_ids_trunc = sender_input_ids[:sender_alloc]
            sender_values_trunc = sender_values[:sender_alloc]
            receiver_input_ids_trunc = receiver_input_ids[:receiver_alloc]
            receiver_values_trunc = receiver_values[:receiver_alloc]

            # Check if truncated sequences are empty
            if not sender_input_ids_trunc or not receiver_input_ids_trunc:
                return None

            # Build cell pair sequence
            pair_input_ids = [cls_token]
            pair_input_ids.extend(sender_input_ids_trunc)
            pair_input_ids.append(sep_token)
            pair_input_ids.extend(receiver_input_ids_trunc)
            pair_input_ids.append(sep_token)

            # Build value sequence
            pair_values = [0.0]
            pair_values.extend(sender_values_trunc)
            pair_values.append(0.0)
            pair_values.extend(receiver_values_trunc)
            pair_values.append(0.0)

            # Pad to fixed length
            current_length = len(pair_input_ids)
            if current_length < max_sequence_length:
                pad_length = max_sequence_length - current_length
                pair_input_ids.extend([pad_token] * pad_length)
                pair_values.extend([0.0] * pad_length)
            else:
                pair_input_ids = pair_input_ids[:max_sequence_length]
                pair_values = pair_values[:max_sequence_length]

            # Final verification
            if len(pair_input_ids) != max_sequence_length or len(pair_values) != max_sequence_length:
                return None

            sequence_length = len(sender_input_ids_trunc) + len(receiver_input_ids_trunc) + 3

            return {
                'input_ids': pair_input_ids,
                'values': pair_values,
                'length': sequence_length
            }

        except Exception as e:
            logger.error(f"Failed to build cell pair sequence: {str(e)}")
            return None

    def predict_interaction_scores(self, dataset, cell_types, max_sequence_length=2048, batch_size=8):
        """
        Predict cell interaction scores
        """
        logger.info("Starting to predict cell interaction scores (regression task)...")

        # Generate all possible cell pairs
        cell_pairs = []
        for i, sender in enumerate(cell_types):
            for j, receiver in enumerate(cell_types):
                if sender != receiver:  # Exclude self-interactions
                    cell_pairs.append((sender, receiver))

        logger.info(f"Need to predict {len(cell_pairs)} cell pairs")

        predictions = []

        with torch.no_grad():
            for i in tqdm(range(0, len(cell_pairs), batch_size), desc="Predicting"):
                batch_predictions = self._predict_batch(
                    dataset, cell_pairs[i:i + batch_size], max_sequence_length
                )
                predictions.extend(batch_predictions)

        logger.info(f"Completed {len(predictions)} predictions")
        return predictions

    def _predict_batch(self, dataset, cell_pairs, max_sequence_length):
        """Batch prediction"""
        batch_predictions = []

        for sender, receiver in cell_pairs:
            try:
                # Build sequence
                sequence = self._build_cell_pair_sequence(dataset, sender, receiver, max_sequence_length)
                if sequence is None:
                    continue

                # Prepare input
                input_ids = torch.tensor([sequence['input_ids']]).long().to(self.device)
                values = torch.tensor([sequence['values']]).float().to(self.device)

                # Predict - regression task directly outputs score
                outputs = self.model(input_ids=input_ids, values=values)
                score = outputs.logits.item()  # Regression task directly takes scalar value

                batch_predictions.append({
                    'sender': sender,
                    'receiver': receiver,
                    'predicted_score': score,
                    'confidence': 1.0  # Regression task has no confidence, set to 1.0
                })
            except Exception as e:
                logger.warning(f"Failed to predict cell pair {sender}-{receiver}: {str(e)}")
                continue

        return batch_predictions

    def create_interaction_matrix(self, predictions, cell_types):
        """
        Create cell interaction scoring matrix
        """
        logger.info("Creating cell interaction scoring matrix...")

        # Create empty DataFrame
        interaction_matrix = pd.DataFrame(
            0.0, index=cell_types, columns=cell_types
        )

        # Fill prediction results
        for pred in predictions:
            sender = pred['sender']
            receiver = pred['receiver']
            score = pred['predicted_score']
            if sender in cell_types and receiver in cell_types:
                interaction_matrix.loc[sender, receiver] = score

        return interaction_matrix

    def save_results(self, interaction_matrix, predictions, output_dir):
        """Save prediction results"""
        os.makedirs(output_dir, exist_ok=True)

        # Save interaction matrix
        matrix_path = os.path.join(output_dir, 'interaction_score_matrix.csv')
        interaction_matrix.to_csv(matrix_path)
        logger.info(f"Interaction scoring matrix saved: {matrix_path}")

        # Save detailed prediction results
        predictions_df = pd.DataFrame(predictions)
        predictions_path = os.path.join(output_dir, 'detailed_predictions.csv')
        predictions_df.to_csv(predictions_path, index=False)
        logger.info(f"Detailed prediction results saved: {predictions_path}")

        # Generate statistical analysis
        stats = self._generate_statistics(interaction_matrix, predictions)
        stats_path = os.path.join(output_dir, 'statistical_analysis.json')
        with open(stats_path, 'w') as f:
            json.dump(stats, f, indent=2)
        logger.info(f"Statistical analysis saved: {stats_path}")

        # Generate visualizations
        self._generate_visualizations(interaction_matrix, output_dir)
        logger.info("Visualization results generated")

    def _generate_statistics(self, interaction_matrix, predictions):
        """Generate statistical analysis"""
        scores = [pred['predicted_score'] for pred in predictions]

        stats = {
            "basic_statistics": {
                "total_cell_types": interaction_matrix.shape[0],
                "total_predicted_interactions": len(predictions),
                "score_range": {
                    "min": float(min(scores)),
                    "max": float(max(scores)),
                    "mean": float(np.mean(scores)),
                    "median": float(np.median(scores)),
                    "std": float(np.std(scores))
                },
                "high_score_interactions": len([s for s in scores if s > np.median(scores)]),
                "low_score_interactions": len([s for s in scores if s <= np.median(scores)])
            },
            "top_interactions": sorted(predictions, key=lambda x: x['predicted_score'], reverse=True)[:20]
        }
        if not predictions:  # Check if predictions is empty
            logger.warning("Prediction results are empty, returning default statistical information")
            return {
                "basic_statistics": {
                    "total_cell_types": interaction_matrix.shape[0],
                    "total_predicted_interactions": 0,
                    "score_range": {
                        "min": 0.0,
                        "max": 0.0,
                        "mean": 0.0,
                        "median": 0.0,
                        "std": 0.0
                    },
                    "high_score_interactions": 0,
                    "low_score_interactions": 0
                },
                "top_interactions": []
            }

        return stats

    def _generate_visualizations(self, interaction_matrix, output_dir):
        """Generate visualization results"""
        try:
            plt.style.use('default')

            # 1. Interaction scoring matrix heatmap
            plt.figure(figsize=(12, 10))
            sns.heatmap(interaction_matrix, annot=True, fmt='.3f', cmap='Reds',
                        cbar_kws={'label': 'Interaction Score'})
            plt.title('Cell Interaction Score Matrix (Regression)')
            plt.tight_layout()
            heatmap_path = os.path.join(output_dir, 'interaction_score_heatmap.png')
            plt.savefig(heatmap_path, dpi=300, bbox_inches='tight')
            plt.close()

            # 2. Score distribution histogram
            plt.figure(figsize=(10, 6))
            scores = interaction_matrix.values.flatten()
            scores = scores[scores > 0]  # Only show positive scores
            plt.hist(scores, bins=50, alpha=0.7, color='skyblue', edgecolor='black')
            plt.xlabel('Interaction Score')
            plt.ylabel('Frequency')
            plt.title('Distribution of Interaction Scores (Regression)')
            plt.grid(True, alpha=0.3)
            hist_path = os.path.join(output_dir, 'score_distribution.png')
            plt.savefig(hist_path, dpi=300, bbox_inches='tight')
            plt.close()

            logger.info("Visualization results generated")
        except Exception as e:
            logger.warning(f"Failed to generate visualization results: {str(e)}")


# Main function
def main():
    """Main execution function"""
    # Configuration parameters - regression task version
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
        'warmup_steps': 100,
        'weight_decay': 0.001,
        'test_size': 0.2,
        'validation_size': 0.1,
        'balance_dataset': False,
        'metric_for_best_model': 'rmse'  # Regression task uses RMSE as best model metric
    }

    # Phase 1: Fine-tune regression model
    logger.info("=== Phase 1: Cell-Cell Interaction Regression Model Fine-Tuning ===")
    trainer, test_predictions, val_predictions = fine_tune_regression_model(config)

    # Print evaluation results
    if trainer is not None and test_predictions is not None:
        logger.info("=== Regression Model Evaluation Results ===")
        for metric, value in test_predictions.metrics.items():
            logger.info(f"{metric}: {value:.4f}")

    # Phase 2: Use fine-tuned regression model for prediction
    if trainer is not None:
        logger.info("=== Phase 2: Cell-Cell Interaction Scoring Prediction ===")

        # Create predictor and perform prediction
        predictor = RegressionCellInteractionPredictor(
            model_path=config['output_dir'],
            token_dict_path=config['token_dict_path']
        )

        # Load dataset to get cell type information
        dataset = load_from_disk(config['dataset_path'])
        if 'cell_type' in dataset.column_names:
            cell_types = sorted(list(set(dataset['cell_type'])))
        else:
            cell_types = [f"Cell_{i}" for i in range(len(dataset))]

        # Perform prediction
        predictions = predictor.predict_interaction_scores(
            dataset=dataset,
            cell_types=cell_types,
            max_sequence_length=config['max_sequence_length'],
            batch_size=config['batch_size']
        )

        # Create interaction scoring matrix
        interaction_matrix = predictor.create_interaction_matrix(predictions, cell_types)

        # Save results
        output_dir = os.path.join(config['output_dir'], 'predictions')
        predictor.save_results(interaction_matrix, predictions, output_dir)

        logger.info("=== Regression Prediction Complete ===")
        logger.info(f"Results saved in: {output_dir}")

    logger.info("=== Processing Complete ===")


if __name__ == "__main__":
    main()
