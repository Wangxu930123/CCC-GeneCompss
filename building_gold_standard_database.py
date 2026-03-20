#!/usr/bin/env python3
"""
Build Gold Standard Database for Cell-Cell Interactions

This script integrates CellChat and CellPhoneDB results to create a gold standard
database for cell-cell interactions.
"""

import os
import pandas as pd
import numpy as np
import matplotlib.pyplot as plt
import seaborn as sns
import networkx as nx
from sklearn.preprocessing import MinMaxScaler
import logging

# Configure logging
logging.basicConfig(
    level=logging.INFO,
    format='%(asctime)s - %(levelname)s - %(message)s',
    datefmt='%Y-%m-%d %H:%M:%S'
)
logger = logging.getLogger(__name__)


def load_cellchat_results(cellchat_dir):
    """
    Load CellChat results
    
    Args:
        cellchat_dir: Directory containing CellChat output files
        
    Returns:
        Dictionary containing CellChat results
    """
    logger.info("Loading CellChat results...")
    results = {}

    # 1. Load interaction strength matrix
    interaction_matrix_path = os.path.join(cellchat_dir, "cell_interaction_strength_matrix.csv")
    if not os.path.exists(interaction_matrix_path):
        alt_names = ['cell_interaction_matrix.csv', 'interaction_matrix.csv', 'interaction_strength.csv']
        for alt in alt_names:
            alt_path = os.path.join(cellchat_dir, alt)
            if os.path.exists(alt_path):
                interaction_matrix_path = alt_path
                break
        else:
            raise FileNotFoundError(f"CellChat interaction matrix file not found: {interaction_matrix_path}")

    logger.info(f"Loading interaction matrix: {interaction_matrix_path}")
    cellchat_matrix = pd.read_csv(interaction_matrix_path)

    # Check column names
    if 'Sender' not in cellchat_matrix.columns:
        sender_cols = [col for col in cellchat_matrix.columns if 'sender' in col.lower() or 'source' in col.lower()]
        if sender_cols:
            cellchat_matrix.rename(columns={sender_cols[0]: 'Sender'}, inplace=True)
        else:
            cellchat_matrix.rename(columns={cellchat_matrix.columns[0]: 'Sender'}, inplace=True)

    # Convert to long format
    receiver_cols = [col for col in cellchat_matrix.columns if col != 'Sender']
    cellchat_matrix = cellchat_matrix.melt(
        id_vars='Sender',
        value_vars=receiver_cols,
        var_name='Receiver',
        value_name='CellChat_Score'
    )

    # Ensure scores are numeric
    cellchat_matrix['CellChat_Score'] = pd.to_numeric(
        cellchat_matrix['CellChat_Score'], errors='coerce'
    ).fillna(0)

    results['matrix'] = cellchat_matrix

    # 2. Load communication results
    communication_path = os.path.join(cellchat_dir, "cellchat_communication.csv")
    if not os.path.exists(communication_path):
        logger.warning(f"CellChat communication file not found: {communication_path}")
        results['communication'] = pd.DataFrame()
        results['pathways'] = pd.DataFrame()
        return results

    logger.info(f"Loading communication file: {communication_path}")
    cellchat_comm = pd.read_csv(communication_path)

    # Rename columns
    col_mapping = {
        'source': ['source', 'sender', 'cell_source', 'from'],
        'target': ['target', 'receiver', 'cell_target', 'to'],
        'interaction_name': ['interaction_name', 'interaction', 'pathway', 'ligand_receptor_pair'],
        'prob': ['prob', 'probability', 'pval', 'p_value', 'score', 'strength']
    }

    for standard_name, alt_names in col_mapping.items():
        for alt in alt_names:
            if alt in cellchat_comm.columns:
                cellchat_comm.rename(columns={alt: standard_name}, inplace=True)
                break

    # Ensure required columns exist
    required_columns = ['source', 'target']
    if not all(col in cellchat_comm.columns for col in required_columns):
        logger.warning(f"Missing required columns in communication file: {[c for c in required_columns if c not in cellchat_comm.columns]}")
        cellchat_comm = pd.DataFrame(columns=['source', 'target', 'interaction_name', 'prob'])
    else:
        if 'interaction_name' not in cellchat_comm.columns:
            cellchat_comm['interaction_name'] = 'Unknown'
        if 'prob' not in cellchat_comm.columns:
            cellchat_comm['prob'] = 0.0

        cellchat_comm['prob'] = pd.to_numeric(cellchat_comm['prob'], errors='coerce').fillna(0)

    cellchat_comm = cellchat_comm[['source', 'target', 'interaction_name', 'prob']]
    cellchat_comm.columns = ['Sender', 'Receiver', 'Interaction', 'CellChat_Prob']
    results['communication'] = cellchat_comm

    return results


def load_cellphonedb_results(cpdb_dir):
    """
    Load CellPhoneDB results
    
    Args:
        cpdb_dir: Directory containing CellPhoneDB output files
        
    Returns:
        DataFrame with aggregated CellPhoneDB results
    """
    logger.info("Loading CellPhoneDB results...")
    
    # Find significant_means file
    sig_means_files = [
        'significant_means.txt', 'significant_means.csv', 'significant_means.tsv',
        'significant_means.xlsx', 'deconvoluted.txt', 'means.txt'
    ]

    sig_means_path = None
    for file_name in sig_means_files:
        test_path = os.path.join(cpdb_dir, file_name)
        if os.path.exists(test_path):
            sig_means_path = test_path
            break

    if not sig_means_path:
        raise FileNotFoundError(f"CellPhoneDB results file not found in directory: {cpdb_dir}")

    logger.info(f"Loading CellPhoneDB file: {sig_means_path}")

    # Try different reading methods
    try:
        if sig_means_path.endswith('.txt') or sig_means_path.endswith('.tsv'):
            cpdb_sig = pd.read_csv(sig_means_path, sep='\t')
        elif sig_means_path.endswith('.csv'):
            cpdb_sig = pd.read_csv(sig_means_path)
        elif sig_means_path.endswith('.xlsx'):
            cpdb_sig = pd.read_excel(sig_means_path)
        else:
            # Auto-detect delimiter
            with open(sig_means_path, 'r') as f:
                first_line = f.readline()

            if '\t' in first_line:
                cpdb_sig = pd.read_csv(sig_means_path, sep='\t')
            elif ',' in first_line:
                cpdb_sig = pd.read_csv(sig_means_path)
            else:
                cpdb_sig = pd.read_csv(sig_means_path, sep=None, engine='python')
    except Exception as e:
        logger.error(f"Failed to parse file: {str(e)}")
        raise

    # Process CellPhoneDB data
    non_pair_cols = ['id_cp_interaction', 'interacting_pair', 'gene_a', 'gene_b', 'partner_a', 'partner_b']
    cell_pair_cols = [col for col in cpdb_sig.columns if col not in non_pair_cols]

    if not cell_pair_cols:
        cell_pair_cols = [col for col in cpdb_sig.columns if '|' in col]

        if not cell_pair_cols:
            numeric_cols = cpdb_sig.select_dtypes(include=np.number).columns
            cell_pair_cols = [col for col in numeric_cols if col not in non_pair_cols]

    if not cell_pair_cols:
        raise ValueError("Cannot identify cell pair columns")

    logger.info(f"Identified {len(cell_pair_cols)} cell pair columns")

    # Convert to long format
    cpdb_long = pd.melt(
        cpdb_sig,
        id_vars=['gene_a', 'gene_b'] if 'gene_a' in cpdb_sig.columns else non_pair_cols,
        value_vars=cell_pair_cols,
        var_name='CellPair',
        value_name='CPDB_Score'
    )

    cpdb_long['CPDB_Score'] = pd.to_numeric(cpdb_long['CPDB_Score'], errors='coerce')
    cpdb_long = cpdb_long.dropna(subset=['CPDB_Score'])
    cpdb_long = cpdb_long[cpdb_long['CPDB_Score'] > 0]

    # Split cell pairs
    try:
        cpdb_long[['Sender', 'Receiver']] = cpdb_long['CellPair'].str.split(r'\|', expand=True)
    except:
        try:
            cpdb_long[['Sender', 'Receiver']] = cpdb_long['CellPair'].str.split(r'[|;:]', expand=True)
        except:
            cpdb_long['Sender'] = cpdb_long['CellPair']
            cpdb_long['Receiver'] = 'Unknown'

    # Aggregate results
    cpdb_agg = cpdb_long.groupby(['Sender', 'Receiver']).agg(
        CPDB_Mean=('CPDB_Score', 'mean'),
        CPDB_Max=('CPDB_Score', 'max'),
        Num_Interactions=('CPDB_Score', 'count')
    ).reset_index()

    return cpdb_agg


def integrate_results(cellchat_results, cpdb_results):
    """
    Integrate CellChat and CellPhoneDB results
    
    Args:
        cellchat_results: Dictionary with CellChat results
        cpdb_results: DataFrame with CellPhoneDB results
        
    Returns:
        DataFrame with integrated results
    """
    logger.info("Integrating results...")
    
    if 'matrix' in cellchat_results and not cellchat_results['matrix'].empty:
        combined = pd.merge(
            cellchat_results['matrix'],
            cpdb_results,
            on=['Sender', 'Receiver'],
            how='outer'
        )
    else:
        combined = cpdb_results.copy()
        combined['CellChat_Score'] = 0

    # Fill missing values
    combined['CellChat_Score'] = combined['CellChat_Score'].fillna(0)
    combined['CPDB_Mean'] = combined['CPDB_Mean'].fillna(0)
    combined['CPDB_Max'] = combined['CPDB_Max'].fillna(0)
    combined['Num_Interactions'] = combined['Num_Interactions'].fillna(0)

    # Normalize scores using MinMaxScaler to [0, 1] range
    scaler = MinMaxScaler()

    if combined['CellChat_Score'].nunique() > 1:
        combined['Norm_CellChat'] = scaler.fit_transform(combined[['CellChat_Score']])
    else:
        combined['Norm_CellChat'] = combined['CellChat_Score']

    if combined['CPDB_Mean'].nunique() > 1:
        combined['Norm_CPDB_Mean'] = scaler.fit_transform(combined[['CPDB_Mean']])
    else:
        combined['Norm_CPDB_Mean'] = combined['CPDB_Mean']

    if combined['CPDB_Max'].nunique() > 1:
        combined['Norm_CPDB_Max'] = scaler.fit_transform(combined[['CPDB_Max']])
    else:
        combined['Norm_CPDB_Max'] = combined['CPDB_Max']

    # Calculate consensus score
    combined['Consensus_Score'] = (
        combined['Norm_CellChat'] +
        combined['Norm_CPDB_Mean'] +
        combined['Norm_CPDB_Max']
    ) / 3
    combined['Consensus_Score'] = combined['Consensus_Score'].fillna(0)

    # Add CellChat communication probability if available
    if 'communication' in cellchat_results and not cellchat_results['communication'].empty:
        cellchat_comm_summary = cellchat_results['communication'].groupby(['Sender', 'Receiver']).agg(
            CellChat_MeanProb=('CellChat_Prob', 'mean'),
            CellChat_MaxProb=('CellChat_Prob', 'max')
        ).reset_index()

        combined = pd.merge(
            combined,
            cellchat_comm_summary,
            on=['Sender', 'Receiver'],
            how='left'
        )

        combined['CellChat_MeanProb'] = combined['CellChat_MeanProb'].fillna(0)
        combined['CellChat_MaxProb'] = combined['CellChat_MaxProb'].fillna(0)
    else:
        combined['CellChat_MeanProb'] = 0
        combined['CellChat_MaxProb'] = 0

    return combined


def generate_full_interaction_matrix(integrated_data):
    """
    Generate complete cell-cell interaction matrix
    
    Args:
        integrated_data: DataFrame with integrated results
        
    Returns:
        DataFrame with complete interaction matrix
    """
    logger.info("Generating complete cell-cell interaction matrix...")

    all_cell_types = sorted(set(integrated_data['Sender'].unique()) | set(integrated_data['Receiver'].unique()))

    # Create all possible cell pair combinations
    all_pairs = pd.DataFrame(
        [(sender, receiver) for sender in all_cell_types for receiver in all_cell_types if sender != receiver],
        columns=['Sender', 'Receiver']
    )

    # Merge with existing data
    full_matrix = pd.merge(
        all_pairs,
        integrated_data,
        on=['Sender', 'Receiver'],
        how='left'
    )

    # Fill missing values
    full_matrix['Consensus_Score'] = full_matrix['Consensus_Score'].fillna(0)
    full_matrix['CellChat_Score'] = full_matrix['CellChat_Score'].fillna(0)
    full_matrix['CPDB_Mean'] = full_matrix['CPDB_Mean'].fillna(0)
    full_matrix['CPDB_Max'] = full_matrix['CPDB_Max'].fillna(0)
    full_matrix['Num_Interactions'] = full_matrix['Num_Interactions'].fillna(0)
    full_matrix['CellChat_MeanProb'] = full_matrix['CellChat_MeanProb'].fillna(0)
    full_matrix['CellChat_MaxProb'] = full_matrix['CellChat_MaxProb'].fillna(0)

    # Add pair ID
    full_matrix['Pair_ID'] = full_matrix['Sender'] + '_' + full_matrix['Receiver']

    return full_matrix


def assign_gold_standard_labels(full_matrix, threshold_method="quantile", threshold_value=0.7):
    """
    Assign gold standard labels to all cell-cell interactions
    
    Args:
        full_matrix: DataFrame with complete interaction matrix
        threshold_method: "quantile" or "absolute"
        threshold_value: Threshold value
        
    Returns:
        DataFrame with gold standard labels
    """
    logger.info("Assigning gold standard labels to all cell-cell interactions...")

    if full_matrix.empty:
        logger.warning("Interaction matrix is empty, cannot assign labels")
        return full_matrix

    labeled_matrix = full_matrix.copy()

    if threshold_method == "quantile":
        threshold = labeled_matrix['Consensus_Score'].quantile(threshold_value)
        logger.info(f"Using quantile threshold: {threshold_value} -> Consensus score threshold: {threshold:.4f}")
    elif threshold_method == "absolute":
        threshold = threshold_value
        logger.info(f"Using absolute threshold: {threshold_value}")
    else:
        raise ValueError("threshold_method must be 'quantile' or 'absolute'")

    # Assign gold standard labels
    # 1: High confidence positive samples (consensus score >= threshold)
    # 0: Negative samples (consensus score < threshold)
    labeled_matrix['Gold_Standard_Label'] = (labeled_matrix['Consensus_Score'] >= threshold).astype(int)

    # Calculate confidence level
    if not labeled_matrix.empty:
        if threshold_method == "quantile":
            high_threshold = threshold
            medium_threshold = labeled_matrix['Consensus_Score'].quantile(threshold_value * 0.7)
            low_threshold = labeled_matrix['Consensus_Score'].quantile(threshold_value * 0.3)
            
            bins = [-np.inf, low_threshold, medium_threshold, high_threshold, np.inf]
            labels = ['Very Low', 'Low', 'Medium', 'High']
        else:
            high_threshold = threshold
            medium_threshold = threshold * 0.7
            low_threshold = threshold * 0.3
            
            bins = [-np.inf, low_threshold, medium_threshold, high_threshold, np.inf]
            labels = ['Very Low', 'Low', 'Medium', 'High']

        labeled_matrix['Confidence_Level'] = pd.cut(
            labeled_matrix['Consensus_Score'],
            bins=bins,
            labels=labels,
            include_lowest=True
        )

        labeled_matrix['Confidence_Level'] = labeled_matrix['Confidence_Level'].astype(str)
        labeled_matrix.loc[labeled_matrix['Gold_Standard_Label'] == 1, 'Confidence_Level'] = 'Gold Standard'
        labeled_matrix['Confidence_Level'] = labeled_matrix['Confidence_Level'].fillna('Very Low')
    else:
        labeled_matrix['Confidence_Level'] = 'Unknown'

    # Statistics
    num_total = len(labeled_matrix)
    num_positive = labeled_matrix['Gold_Standard_Label'].sum()
    num_negative = num_total - num_positive

    logger.info(f"Total interactions: {num_total}")
    logger.info(f"Gold standard positive samples: {num_positive} ({num_positive / num_total * 100:.1f}%)")
    logger.info(f"Negative samples: {num_negative} ({num_negative / num_total * 100:.1f}%)")

    return labeled_matrix


def visualize_all_interactions(labeled_matrix, output_dir):
    """
    Visualize all cell-cell interactions
    
    Args:
        labeled_matrix: DataFrame with labeled interactions
        output_dir: Directory to save visualizations
    """
    logger.info("Visualizing all cell-cell interactions...")
    os.makedirs(output_dir, exist_ok=True)

    if labeled_matrix.empty:
        logger.warning("Labeled matrix is empty, skipping visualization")
        return

    # 1. Complete interaction matrix heatmap
    heatmap_data = labeled_matrix.pivot_table(
        index='Sender',
        columns='Receiver',
        values='Consensus_Score',
        fill_value=0
    )

    heatmap_data = heatmap_data.reindex(
        index=sorted(heatmap_data.index),
        columns=sorted(heatmap_data.columns)
    )

    plt.figure(figsize=(16, 14))
    sns.heatmap(heatmap_data, cmap='Reds', annot=True, fmt=".3f", linewidths=.5,
                cbar_kws={'label': 'Consensus Score'})
    plt.title('Complete Cell-Cell Interaction Matrix (All Pairs)')
    plt.tight_layout()
    plt.savefig(os.path.join(output_dir, 'complete_interaction_matrix.png'), dpi=300, bbox_inches='tight')
    plt.close()

    # 2. Gold standard label heatmap
    label_heatmap = labeled_matrix.pivot_table(
        index='Sender',
        columns='Receiver',
        values='Gold_Standard_Label',
        fill_value=0
    )

    label_heatmap = label_heatmap.reindex(
        index=sorted(label_heatmap.index),
        columns=sorted(label_heatmap.columns)
    )

    plt.figure(figsize=(16, 14))
    sns.heatmap(label_heatmap, cmap=['lightgray', 'red'], annot=True, fmt="d",
                linewidths=.5, cbar_kws={'label': 'Gold Standard Label (0/1)'})
    plt.title('Gold Standard Labels for All Cell-Cell Interactions')
    plt.tight_layout()
    plt.savefig(os.path.join(output_dir, 'gold_standard_labels.png'), dpi=300, bbox_inches='tight')
    plt.close()

    logger.info("Visualizations saved successfully")


def create_machine_learning_dataset(labeled_matrix, include_features=True):
    """
    Create dataset for machine learning
    
    Args:
        labeled_matrix: DataFrame with labeled interactions
        include_features: Whether to include feature columns
        
    Returns:
        DataFrame with machine learning dataset
    """
    logger.info("Creating machine learning dataset...")

    if labeled_matrix.empty:
        logger.warning("Labeled matrix is empty, cannot create dataset")
        return pd.DataFrame()

    ml_dataset = labeled_matrix[['Sender', 'Receiver', 'Pair_ID', 'Gold_Standard_Label', 'Confidence_Level']].copy()

    if include_features:
        feature_columns = [
            'Consensus_Score', 'CellChat_Score', 'CPDB_Mean', 'CPDB_Max',
            'Num_Interactions', 'CellChat_MeanProb', 'CellChat_MaxProb'
        ]

        available_features = [col for col in feature_columns if col in labeled_matrix.columns]
        for feature in available_features:
            ml_dataset[feature] = labeled_matrix[feature]

        logger.info(f"Included features: {', '.join(available_features)}")

    # Add cell type encoding
    all_cell_types = sorted(set(labeled_matrix['Sender'].unique()) | set(labeled_matrix['Receiver'].unique()))
    cell_type_map = {cell_type: idx for idx, cell_type in enumerate(all_cell_types)}

    ml_dataset['Sender_Encoded'] = ml_dataset['Sender'].map(cell_type_map)
    ml_dataset['Receiver_Encoded'] = ml_dataset['Receiver'].map(cell_type_map)

    logger.info(f"Total cell types: {len(all_cell_types)}")
    logger.info(f"Machine learning dataset shape: {ml_dataset.shape}")

    return ml_dataset


def create_complete_gold_standard_dataset(cellchat_dir, cpdb_dir, output_dir,
                                      threshold_method="quantile", threshold_value=0.7):
    """
    Create complete gold standard dataset for cell-cell interactions
    
    Args:
        cellchat_dir: Directory containing CellChat results
        cpdb_dir: Directory containing CellPhoneDB results
        output_dir: Directory to save output files
        threshold_method: "quantile" or "absolute"
        threshold_value: Threshold value
        
    Returns:
        Dictionary with results
    """
    # Create output directory
    os.makedirs(output_dir, exist_ok=True)

    try:
        cellchat_results = load_cellchat_results(cellchat_dir)
    except Exception as e:
        logger.error(f"Error loading CellChat results: {str(e)}")
        return None

    try:
        cpdb_results = load_cellphonedb_results(cpdb_dir)
    except Exception as e:
        logger.error(f"Error loading CellPhoneDB results: {str(e)}")
        return None

    # Integrate results
    integrated_data = integrate_results(cellchat_results, cpdb_results)

    # Generate complete interaction matrix
    full_matrix = generate_full_interaction_matrix(integrated_data)

    # Assign gold standard labels
    labeled_matrix = assign_gold_standard_labels(
        full_matrix,
        threshold_method=threshold_method,
        threshold_value=threshold_value
    )

    # Create machine learning dataset
    ml_dataset = create_machine_learning_dataset(labeled_matrix, include_features=True)

    # Save all results
    full_output_path = os.path.join(output_dir, 'complete_labeled_interactions.csv')
    labeled_matrix.to_csv(full_output_path, index=False)
    logger.info(f"Complete labeled matrix saved: {full_output_path}")

    ml_output_path = os.path.join(output_dir, 'machine_learning_dataset.csv')
    ml_dataset.to_csv(ml_output_path, index=False)
    logger.info(f"Machine learning dataset saved: {ml_output_path}")

    gold_standard_only = labeled_matrix[labeled_matrix['Gold_Standard_Label'] == 1]
    gold_output_path = os.path.join(output_dir, 'gold_standard_interactions.csv')
    gold_standard_only.to_csv(gold_output_path, index=False)
    logger.info(f"Gold standard interactions saved: {gold_output_path}")

    # Visualize results
    visualize_all_interactions(labeled_matrix, output_dir)

    # Save dataset statistics
    stats = {
        'total_interactions': len(labeled_matrix),
        'gold_standard_positive': len(gold_standard_only),
        'gold_standard_negative': len(labeled_matrix) - len(gold_standard_only),
        'positive_percentage': len(gold_standard_only) / len(labeled_matrix) * 100,
        'threshold_method': threshold_method,
        'threshold_value': threshold_value,
        'num_cell_types': len(set(labeled_matrix['Sender'].unique()) | set(labeled_matrix['Receiver'].unique()))
    }

    stats_df = pd.DataFrame([stats])
    stats_path = os.path.join(output_dir, 'dataset_statistics.csv')
    stats_df.to_csv(stats_path, index=False)
    logger.info(f"Dataset statistics saved: {stats_path}")

    logger.info(f"\n✅ Gold standard dataset created successfully!")
    logger.info(f"Total interactions: {len(labeled_matrix)}")
    logger.info(f"Gold standard positive samples: {len(gold_standard_only)}")
    logger.info(f"Positive percentage: {stats['positive_percentage']:.1f}%")
    logger.info(f"Cell types: {stats['num_cell_types']}")
    logger.info(f"Threshold method: {threshold_method}, Threshold: {threshold_value}")

    return {
        'labeled_matrix': labeled_matrix,
        'ml_dataset': ml_dataset,
        'gold_standard_only': gold_standard_only,
        'statistics': stats
    }


if __name__ == "__main__":
    import argparse
    
    parser = argparse.ArgumentParser(description='Build gold standard database for cell-cell interactions')
    parser.add_argument('--cellchat_dir', type=str, required=True,
                        help='Directory containing CellChat results')
    parser.add_argument('--cpdb_dir', type=str, required=True,
                        help='Directory containing CellPhoneDB results')
    parser.add_argument('--output_dir', type=str, required=True,
                        help='Directory to save output files')
    parser.add_argument('--threshold_method', type=str, default='quantile',
                        choices=['quantile', 'absolute'],
                        help='Threshold method: quantile or absolute')
    parser.add_argument('--threshold_value', type=float, default=0.7,
                        help='Threshold value (0-1 for quantile, absolute score for absolute)')
    
    args = parser.parse_args()
    
    create_complete_gold_standard_dataset(
        cellchat_dir=args.cellchat_dir,
        cpdb_dir=args.cpdb_dir,
        output_dir=args.output_dir,
        threshold_method=args.threshold_method,
        threshold_value=args.threshold_value
    )
