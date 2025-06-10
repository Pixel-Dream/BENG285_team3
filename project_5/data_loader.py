import scanpy as sc
import pandas as pd
import numpy as np
from typing import Optional, Literal

def load_and_preprocess_data(file_path: str) -> sc.AnnData:
    """Load h5ad file and filter out samples with NaN survival labels."""
    dataset = sc.read_h5ad(file_path)
    dataset = dataset[dataset.obs['survival_3yr_label'].notna()].copy()
    dataset.obs['survival_3yr_label'] = dataset.obs['survival_3yr_label'].astype(int)
    return dataset

def construct_feature_matrix(
    adata: sc.AnnData,
    expression_type: str = 'normalized',
    mutation_data_source: Optional[Literal['anndata', 'csv']] = 'anndata',
    mutation_data_identifier: Optional[str] = 'jiaming_data_full',
) -> pd.DataFrame:
    """
    Construct a feature matrix from AnnData object combining gene expression and mutation data.
    """
    sample_ids = adata.obs_names.tolist()
    
    # 1. Get gene expression data
    if expression_type == 'normalized':
        expr_data = adata.X
    elif expression_type == 'raw':
        if 'counts' not in adata.layers:
            raise ValueError("No 'counts' layer found in adata.layers")
        expr_data = adata.layers['counts']
    else:
        raise ValueError(f"expression_type must be 'normalized' or 'raw', got {expression_type}")
    
    if hasattr(expr_data, 'toarray'):
        expr_data = expr_data.toarray()

    gene_names = adata.var_names.tolist()
    expr_df = pd.DataFrame(
        expr_data,
        index=sample_ids,
        columns=[f'gene_{gene}' for gene in gene_names]
    )
    
    # 2. Add mutation data if requested
    if mutation_data_source:
        if mutation_data_source == 'anndata':
            if mutation_data_identifier not in adata.uns:
                raise ValueError(f"'{mutation_data_identifier}' not found in adata.uns")
            mut_df = adata.uns[mutation_data_identifier].copy()
        
        elif mutation_data_source == 'csv':
            try:
                mut_df = pd.read_csv(mutation_data_identifier, index_col=0)
            except FileNotFoundError:
                raise FileNotFoundError(f"Mutation CSV file not found at: {mutation_data_identifier}")
        else:
            raise ValueError(f"Unknown mutation_data_source: {mutation_data_source}")

        mut_df.columns = [f'mut_{col}' for col in mut_df.columns]
        
        # Align mutation data with expression data
        feature_matrix = expr_df.join(mut_df, how='left').fillna(0)
    else:
        feature_matrix = expr_df

    return feature_matrix

def get_feature_counts(feature_matrix: pd.DataFrame) -> dict:
    """Count the number of DNA (mutation) and RNA (gene expression) features."""
    feature_names = feature_matrix.columns.tolist()
    
    rna_features = [f for f in feature_names if f.startswith('gene_')]
    dna_features = [f for f in feature_names if f.startswith('mut_')]
    
    return {
        'n_rna_features': len(rna_features),
        'n_dna_features': len(dna_features),
        'total_features': len(feature_names),
        'rna_feature_names': rna_features,
        'dna_feature_names': dna_features,
    } 