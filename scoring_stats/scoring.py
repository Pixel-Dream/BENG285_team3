#!/usr/bin/env python3
# embedding_clustering_evaluator.py
"""
This script evaluates the effectiveness of different embedding methods for clustering data.
It compares how well different clustering algorithms perform on various embeddings and
provides metrics to assess clustering quality.

The script uses true labels from an AnnData object for evaluation.
"""

import os
import pandas as pd
import numpy as np
import matplotlib.pyplot as plt
import seaborn as sns
from sklearn.cluster import KMeans
import hdbscan
from sklearn.metrics import adjusted_rand_score, adjusted_mutual_info_score, silhouette_score
from pathlib import Path
import argparse
import scanpy as sc
import warnings
warnings.filterwarnings('ignore')


def load_embeddings(base_dir):
    """
    Load all embedding CSV files from directories
    
    Args:
        base_dir (str): Base directory containing embedding subdirectories
        
    Returns:
        dict: Dictionary with embedding method names as keys and dataframes as values
    """
    embedding_files = {}
    for root, _, files in os.walk(base_dir):
        for file in files:
            if file == 'embedding_coordinates.csv':  # Exact filename match
                path = Path(root) / file
                method_name = path.parent.name
                df = pd.read_csv(path, index_col=0)
                embedding_files[method_name] = df
    
    print(f"Loaded {len(embedding_files)} embedding methods: {list(embedding_files.keys())}")
    return embedding_files


def extract_labels_from_anndata(anndata_path, label_columns=None):
    """
    Extract true labels from an AnnData object
    
    Args:
        anndata_path (str): Path to AnnData object (.h5ad file)
        label_columns (list): List of column names to extract as labels
        
    Returns:
        dict: Dictionary with label_name:label_series pairs
    """
    print(f"Loading AnnData from {anndata_path}")
    adata = sc.read_h5ad(anndata_path)
    
    # If no label columns specified, use defaults that are likely to be of interest
    if label_columns is None:
        label_columns = [
            'demographic.gender', 
            'smoker_status_grouped',
            'tumor_status'
        ]
    
    # Verify columns exist
    available_columns = adata.obs.columns
    valid_columns = [col for col in label_columns if col in available_columns]
    
    if len(valid_columns) == 0:
        raise ValueError(f"None of the requested columns {label_columns} found in AnnData. "
                         f"Available columns: {list(available_columns)}")
    
    # Extract labels
    labels = {}
    for col in valid_columns:
        # Convert category to string to avoid issues
        if pd.api.types.is_categorical_dtype(adata.obs[col]):
            labels[col] = adata.obs[col].astype(str)
        else:
            labels[col] = adata.obs[col]
            
        # Print label distribution
        print(f"Label distribution for {col}:")
        print(labels[col].value_counts())
        print()
    
    return labels


def perform_clustering(df, n_clusters=None, min_cluster_size=5, min_samples=5, random_state=42):
    """
    Perform multiple clustering algorithms on the given data
    
    Args:
        df (DataFrame): Embedding coordinates dataframe
        n_clusters (dict): Dictionary with clustering method names and number of clusters
        min_cluster_size (int): Minimum cluster size for HDBSCAN
        min_samples (int): Minimum samples parameter for HDBSCAN
        random_state (int): Random state for reproducibility
        
    Returns:
        dict: Dictionary with clustering method names as keys and cluster labels as values
    """
    results = {}
    
    # Get only the numeric columns for clustering
    numeric_cols = df.select_dtypes(include=np.number).columns
    data = df[numeric_cols].values
    
    # Set default n_clusters if not provided
    if n_clusters is None:
        n_clusters = {
            'kmeans_gender': 2,
            'kmeans_smoker': 3,
            'kmeans_tumor': 4
        }
    
    # KMeans with different numbers of clusters for different label types
    for method, n in n_clusters.items():
        kmeans = KMeans(n_clusters=n, random_state=random_state)
        results[method] = kmeans.fit_predict(data)
    
    # HDBSCAN - density-based clustering
    hdb = hdbscan.HDBSCAN(min_cluster_size=min_cluster_size, min_samples=min_samples)
    results['hdbscan'] = hdb.fit_predict(data)
    
    return results


def calculate_metrics(clustering_results, true_labels_dict, df):
    """
    Calculate clustering evaluation metrics
    
    Args:
        clustering_results (dict): Dictionary with clustering method names and cluster labels
        true_labels_dict (dict): Dictionary with label names and label series
        df (DataFrame): Original dataframe with embedding coordinates
        
    Returns:
        DataFrame: Metrics for each clustering method and true label
    """
    metrics = []
    
    # Get only the numeric columns for clustering
    numeric_cols = df.select_dtypes(include=np.number).columns
    data = df[numeric_cols].values
    
    for method, labels in clustering_results.items():
        # Skip if all samples are assigned to noise (-1) or only one cluster
        if len(np.unique(labels[labels != -1])) <= 1:
            continue
            
        # Calculate silhouette score across all data points
        silhouette = np.nan
        valid_indices = labels != -1
        
        if np.sum(valid_indices) > 0 and len(np.unique(labels[valid_indices])) > 1:
            silhouette = silhouette_score(
                data[valid_indices], 
                labels[valid_indices], 
                metric='euclidean'
            )
            
        # Compare to each true label type
        for label_name, true_labels in true_labels_dict.items():
            # --- START Check for matching sample IDs ---
            if not df.index.equals(true_labels.index):
                raise ValueError(
                    f"Sample IDs mismatch between embedding data and true labels for '{label_name}'. "
                    "Ensure both are indexed and ordered identically."
                )
            # --- END Check ---
            
            # Calculate metrics
            ari = adjusted_rand_score(true_labels, labels)
            ami = adjusted_mutual_info_score(true_labels, labels)
            
            metrics.append({
                'method': method,
                'true_label_type': label_name,
                'ari': ari,
                'ami': ami,
                'silhouette': silhouette,
                'n_clusters': len(np.unique(labels[labels != -1])),
                'noise_points': np.sum(labels == -1),
                'noise_percent': np.sum(labels == -1) / len(labels) * 100
            })
    
    return pd.DataFrame(metrics)


def plot_results(results, output_dir='results'):
    """
    Plot evaluation results
    
    Args:
        results (DataFrame): Results from clustering evaluation
        output_dir (str): Directory to save plot results
    """
    os.makedirs(output_dir, exist_ok=True)
    
    # Group results by true label type for separate plots
    for label_type in results['true_label_type'].unique():
        label_results = results[results['true_label_type'] == label_type]
        
        # Create heatmaps for each metric
        for metric in ['ari', 'ami', 'silhouette']:
            pivot_df = label_results.pivot(
                index='embedding_method', 
                columns='method', 
                values=metric
            )
            
            plt.figure(figsize=(12, 8))
            sns.heatmap(
                pivot_df, 
                annot=True, 
                cmap='viridis', 
                fmt='.3f',
                linewidths=.5
            )
            plt.title(f'{metric.upper()} Score by Embedding Method and Clustering Algorithm\nTrue Label: {label_type}')
            plt.tight_layout()
            plt.savefig(os.path.join(output_dir, f'{label_type}_{metric}_heatmap.png'))
            plt.close()
        
        # Bar charts for each metric
        metrics = ['ari', 'ami', 'silhouette']
        for metric in metrics:
            plt.figure(figsize=(14, 8))
            ax = sns.barplot(
                data=label_results, 
                x='embedding_method', 
                y=metric, 
                hue='method'
            )
            plt.title(f'{metric.upper()} by Embedding Method and Clustering Algorithm\nTrue Label: {label_type}')
            plt.xticks(rotation=45, ha='right')
            plt.legend(title='Clustering Method')
            plt.tight_layout()
            plt.savefig(os.path.join(output_dir, f'{label_type}_{metric}_barplot.png'))
            plt.close()
    
    # Create a summary bar chart showing the best performing method for each label type
    plt.figure(figsize=(15, 10))
    summary = results.groupby(['true_label_type', 'embedding_method'])['ari'].mean().reset_index()
    sns.barplot(data=summary, x='true_label_type', y='ari', hue='embedding_method')
    plt.title('Average ARI Score by True Label Type and Embedding Method')
    plt.xlabel('True Label Type')
    plt.ylabel('ARI Score')
    plt.xticks(rotation=45, ha='right')
    plt.tight_layout()
    plt.savefig(os.path.join(output_dir, 'summary_ari_by_label_type.png'))
    plt.close()


def evaluate_embeddings(base_dir, anndata_path, label_columns=None, output_dir='results'):
    """
    Main function to evaluate all embeddings
    
    Args:
        base_dir (str): Base directory containing embedding subdirectories
        anndata_path (str): Path to the AnnData object (.h5ad file)
        label_columns (list): List of column names to extract as labels
        output_dir (str): Directory to save results
        
    Returns:
        DataFrame: Results of the evaluation
    """
    # Create output directory
    os.makedirs(output_dir, exist_ok=True)
    
    # Extract true labels from AnnData
    true_labels_dict = extract_labels_from_anndata(anndata_path, label_columns)
    
    # Load embeddings
    embeddings = load_embeddings(base_dir)
    
    # Results storage
    all_results = []
    
    # Define number of clusters for each clustering task based on label information
    n_clusters = {}
    for label_name, label_series in true_labels_dict.items():
        n_unique = len(label_series.unique())
        method_name = f"kmeans_{label_name.split('.')[-1]}"  # Extract last part of column name
        n_clusters[method_name] = n_unique
    
    # Process each embedding
    for method, df in embeddings.items():
        print(f"Processing embedding method: {method}")
        
        # Perform clustering
        clustering_results = perform_clustering(df, n_clusters=n_clusters)
        
        # Calculate metrics
        metrics = calculate_metrics(clustering_results, true_labels_dict, df)
        metrics['embedding_method'] = method
        
        all_results.append(metrics)
    
    # Combine results
    if all_results:
        final_results = pd.concat(all_results, ignore_index=True)
        
        # Save results to CSV
        results_path = os.path.join(output_dir, 'clustering_evaluation_results.csv')
        final_results.to_csv(results_path, index=False)
        
        # Plot results
        plot_results(final_results, output_dir)
        
        # Print best methods
        print("\nBest methods by true label type:")
        for label_type in final_results['true_label_type'].unique():
            label_df = final_results[final_results['true_label_type'] == label_type]
            best_ari = label_df.loc[label_df['ari'].idxmax()]
            print(f"{label_type}: {best_ari['embedding_method']} with {best_ari['method']} "
                  f"(ARI: {best_ari['ari']:.3f}, AMI: {best_ari['ami']:.3f})")
        
        return final_results
    else:
        print("No results were generated. Check your embedding files and paths.")
        return None


def main():
    """
    Main entry point for the script
    """
    parser = argparse.ArgumentParser(description='Evaluate embeddings for clustering')
    parser.add_argument('--base_dir', type=str, required=True,
                        help='Base directory containing embedding subdirectories')
    parser.add_argument('--anndata', type=str, required=True,
                        help='Path to AnnData object (.h5ad file)')
    parser.add_argument('--output_dir', type=str, default='results',
                        help='Directory to save results (default: results)')
    parser.add_argument('--label_columns', type=str, nargs='+',
                        help='Columns to use as true labels (optional, defaults to gender, smoker_status_grouped, tumor_status)')
    
    args = parser.parse_args()
    
    # Run the evaluation
    results = evaluate_embeddings(
        args.base_dir,
        args.anndata,
        args.label_columns,
        args.output_dir
    )
    
    if results is not None:
        print(f"\nEvaluation complete! Results saved to {args.output_dir}")


if __name__ == "__main__":
    main()