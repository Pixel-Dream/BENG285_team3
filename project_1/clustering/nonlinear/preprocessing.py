import numpy as np
import pandas as pd
import scanpy as sc
import anndata as ad
import json

def preprocess_metadata(adata, group_mapping, inplace=False):
    """
    Preprocess AnnData metadata (obs) based on provided group mappings
    
    Parameters:
    -----------
    adata : AnnData
        AnnData object with original metadata
    group_mapping : dict or str
        Dictionary with grouping information or path to JSON file
    inplace : bool
        Whether to modify the original AnnData object or return a copy
        
    Returns:
    --------
    adata : AnnData
        AnnData object with updated metadata
    """
    # If not inplace, create a copy of the AnnData object
    if not inplace:
        adata = adata.copy()
    
    # Load grouping information if provided as a file path
    if isinstance(group_mapping, str):
        with open(group_mapping, 'r') as f:
            group_mapping = json.load(f)
    
    # For tracking the changes made
    updated_columns = []
    
    # Process each category in the mapping
    for category, fields in group_mapping.items():
        for field, group_values in fields.items():
            # Construct the obs key (category.field)
            obs_key = f"{category}.{field}"
            
            # Check if the key exists in the obs dataframe
            if obs_key not in adata.obs.columns:
                print(f"Warning: {obs_key} not found in AnnData obs")
                continue
            
            # Create a new column name for the grouped data
            grouped_col = f"{obs_key}_grouped"
            
            # Initialize with 'others' for all values
            adata.obs[grouped_col] = 'others'
            
            # Process each group value
            for group_value in group_values:
                # Handle special case of slash-separated values
                if '/' in group_value:
                    # Split by slash and create a list of individual values
                    individual_values = group_value.split('/')
                    
                    # Create mask for values matching any individual value
                    mask = adata.obs[obs_key].isin(individual_values)
                    
                    # Update grouped column where mask is True
                    adata.obs.loc[mask, grouped_col] = group_value
                else:
                    # Simple case: direct value matching
                    mask = adata.obs[obs_key] == group_value
                    adata.obs.loc[mask, grouped_col] = group_value
            
            # Track the updated column
            updated_columns.append(grouped_col)
    
    # Print summary of preprocessing
    print(f"Metadata preprocessing complete. Created {len(updated_columns)} grouped columns:")
    for col in updated_columns:
        value_counts = adata.obs[col].value_counts()
        print(f"  - {col}: {len(value_counts)} groups")
        # Print distribution of values for each grouped column
        for val, count in value_counts.items():
            print(f"    - {val}: {count} observations ({count/len(adata.obs)*100:.1f}%)")
    
    return adata


def preprocess_anndata(adata, n_top_genes=None):
    """
    Preprocess an AnnData object for t-SNE analysis
    
    Parameters:
    -----------
    adata : AnnData
        AnnData object containing expression data and metadata
    n_top_genes : int or None
        Number of most variable genes to select; if None, use all
        
    Returns:
    --------
    expr_data : DataFrame
        Processed expression data with samples as rows
    metadata : DataFrame
        Metadata with samples as index
    """
    print(f"Working with AnnData object: {adata.shape[0]} samples, {adata.shape[1]} genes")
    
    # Create a copy to avoid modifying the original
    adata_copy = adata.copy()
    
    # Extract metadata
    metadata = adata_copy.obs.copy()
    
    # Select highly variable genes if specified
    if n_top_genes is not None and n_top_genes < adata_copy.shape[1]:
        print(f"Selecting {n_top_genes} most variable genes...")
        
        # Check if variable genes have been calculated
        if 'vst.variance.standardized' in adata_copy.var.columns:
            # Use pre-calculated variance
            print("Using pre-calculated gene variance")
            gene_vars = adata_copy.var['vst.variance.standardized']
        else:
            # Calculate highly variable genes
            print("Calculating gene variances")
            sc.pp.highly_variable_genes(adata_copy, n_top_genes=n_top_genes, flavor='seurat')
            if 'highly_variable' in adata_copy.var:
                gene_vars = adata_copy.var['dispersions_norm']
            else:
                # Compute variance manually
                gene_vars = np.var(adata_copy.X, axis=0)
                if isinstance(gene_vars, np.matrix):
                    gene_vars = gene_vars.A1
        
        # Select top genes
        top_gene_idx = np.argsort(gene_vars)[-n_top_genes:]
        adata_copy = adata_copy[:, top_gene_idx]
    
    # Convert to dense if sparse
    # TODO: CHECK HERE
    if isinstance(adata_copy.X, np.matrix) or hasattr(adata_copy.X, 'toarray'):
        expr_matrix = adata_copy.X.toarray()
    else:
        expr_matrix = adata_copy.X
    
    # For expression data
    expr_data = pd.DataFrame(
        expr_matrix, 
        index=adata_copy.obs_names,
        columns=adata_copy.var_names
    )
    
    return expr_data, metadata


def get_metadata_types(metadata, exclude_mostly_null=True, null_threshold=0.8):
    """
    Auto-detect metadata types (categorical vs continuous)
    
    Parameters:
    -----------
    metadata : DataFrame
        Metadata DataFrame
    exclude_mostly_null : bool
        Whether to exclude columns with too many null values
    null_threshold : float
        Threshold for excluding columns (fraction of null values)
        
    Returns:
    --------
    metadata_types : dict
        Dictionary mapping column names to types (categorical/continuous)
    """
    metadata_types = {}
    
    for col in metadata.columns:
        # Check for null values
        null_fraction = metadata[col].isna().mean()
        if exclude_mostly_null and null_fraction > null_threshold:
            continue
        
        # Check data type
        if metadata[col].dtype.name in ['object', 'category', 'bool']:
            metadata_types[col] = 'categorical'
        elif metadata[col].nunique() < 10:
            metadata_types[col] = 'categorical'
        else:
            metadata_types[col] = 'continuous'
    
    return metadata_types