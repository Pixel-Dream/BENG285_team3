import os
import numpy as np
import pandas as pd
import matplotlib.pyplot as plt
import itertools

from preprocessing import preprocess_metadata, preprocess_anndata, get_metadata_types
from visualization import plot_and_save_embedding
from utils import create_parameter_dir, create_method_factory

# Create method factory
get_reduction_method = create_method_factory()

def run_dimensionality_reduction_with_plots(adata, method_name='tsne', params=None,
                                          group_mapping=None, output_prefix='result',
                                          n_top_genes=2000, color_by=None, categorical=None,
                                          output_dir='plots', random_state=42):
    """
    Run dimensionality reduction and generate individual plots for each metadata field
    
    Parameters:
    -----------
    adata : AnnData
        AnnData object with expression data and metadata
    method_name : str
        Name of the dimensionality reduction method ('tsne', 'umap')
    params : dict or None
        Dictionary with method parameters (if None, use defaults)
    group_mapping : dict or str
        Dictionary with grouping information or path to JSON file
    output_prefix : str
        Prefix for output files
    n_top_genes : int
        Number of most variable genes to use
    color_by : list or None
        Metadata columns to color by, if None, plot all, if [], plot one without color
    categorical : list or None
        Boolean flags for each color_by column (if None, auto-detect)
    output_dir : str
        Directory for saving plots
    random_state : int or None
        Random state seed for reproducibility (used by methods like t-SNE, UMAP)
        
    Returns:
    --------
    result : DataFrame
        DataFrame with embedding coordinates
    metadata : DataFrame
        Metadata used for visualization
    adata_grouped : AnnData
        AnnData object with grouped metadata
    plot_paths : list
        Paths to saved plot files
    """
    # First, preprocess the metadata to create grouped columns
    adata_grouped = preprocess_metadata(adata, group_mapping, inplace=False)
    
    # Preprocess AnnData for dimensionality reduction
    expr_data, metadata = preprocess_anndata(
        adata_grouped, 
        n_top_genes=n_top_genes
    )
    
    # Get dimensionality reduction method
    reduction_method = get_reduction_method(method_name)
    
    # Use default parameters if not specified
    if params is None:
        params = reduction_method.get_default_params()
    
    # Ensure random_state is included in params if the method supports it
    params_with_random_state = params.copy()
    params_with_random_state['random_state'] = random_state
    
    # Run dimensionality reduction
    result = reduction_method.fit_transform(expr_data, **params_with_random_state)
    
    # Save results
    method_dir = create_parameter_dir(output_dir, reduction_method.get_method_name(), params)
    result_path = os.path.join(method_dir, f"{output_prefix}_coordinates.csv")
    result.to_csv(result_path)
    print(f"Saved embedding coordinates to: {result_path}")
    
    # Auto-detect metadata types
    metadata_types = get_metadata_types(metadata)
    
    # Auto-select interesting metadata if not specified
    if color_by is None:
        # Find all grouped columns
        color_by = [col for col in metadata.columns if col.endswith('_grouped')]
    
    # Prepare categorical flags if not provided
    if categorical is None:
        categorical = [metadata_types.get(col, 'categorical') == 'categorical' 
                      for col in color_by]
    
    # Generate individual plots
    plot_paths = plot_and_save_embedding(
        result, metadata, 
        method_name=reduction_method.get_method_name(),
        color_by=color_by,
        categorical=categorical,
        params=params,
        output_dir=output_dir
    )
    
    # Return results
    return result, metadata, adata_grouped, plot_paths, result_path


def explore_parameters_with_plots(adata, method_name='tsne', param_grid=None,
                                group_mapping=None, color_by=None, categorical=None,
                                n_top_genes=2000, output_dir='plots', random_state=42):
    """
    Explore parameters for a dimensionality reduction method
    
    Parameters:
    -----------
    adata : AnnData
        AnnData object with expression data and metadata
    method_name : str
        Name of the dimensionality reduction method ('tsne', 'umap')
    param_grid : dict or None
        Dictionary mapping parameter names to lists of values (if None, use default grid)
    group_mapping : dict or str
        Dictionary with grouping information or path to JSON file
    color_by : list or None
        Metadata columns to color by (if None, auto-select interesting grouped ones)
    categorical : list or None
        Boolean flags for each color_by column (if None, auto-detect)
    n_top_genes : int
        Number of most variable genes to use
    output_dir : str
        Directory for saving plots
    random_state : int or None
        Random state seed for reproducibility (used by methods like t-SNE, UMAP)
        
    Returns:
    --------
    results : dict
        Dictionary of results for different parameter combinations
    adata_grouped : AnnData
        AnnData object with grouped metadata
    plot_paths : dict
        Dictionary mapping parameter combinations to lists of plot paths
    hyperparameter_df : DataFrame
        DataFrame with hyperparameters and file paths to saved embeddings
    """
    # First, preprocess the metadata to create grouped columns
    adata_grouped = preprocess_metadata(adata, group_mapping, inplace=False)
    
    # Preprocess data once
    expr_data, metadata = preprocess_anndata(
        adata_grouped, 
        n_top_genes=n_top_genes
    )
    
    # Auto-detect metadata types
    metadata_types = get_metadata_types(metadata)
    
    # Auto-select interesting metadata if not specified
    if color_by is None:
        # Find all grouped columns
        grouped_columns = [col for col in metadata.columns if col.endswith('_grouped')]
        
        # Use grouped columns as default if available
        if grouped_columns:
            color_by = grouped_columns[:5]  # Limit to 5 fields for parameter exploration
        else:
            # Suggested fields from the original
            default_fields = [
                'diagnoses.primary_diagnosis', 
                'diagnoses.ajcc_pathologic_stage',
                'demographic.gender',
                'exposures.tobacco_smoking_status'
            ]
            
            color_by = [field for field in default_fields if field in metadata.columns]
            
            # Add a few more if available
            if len(color_by) < 3:
                additional_fields = list(metadata_types.keys())[:5]
                color_by.extend([f for f in additional_fields if f not in color_by])
                color_by = color_by[:5]  # Limit to 5 fields
    
    # Prepare categorical flags if not provided
    if categorical is None:
        categorical = [metadata_types.get(col, 'categorical') == 'categorical' 
                       for col in color_by]
    
    # Get dimensionality reduction method
    reduction_method = get_reduction_method(method_name)
    
    # Get default parameter grid if not specified
    if param_grid is None:
        param_grid = reduction_method.get_param_grid()
    
    # Create base output directory
    os.makedirs(output_dir, exist_ok=True)
    
    # Initialize results dictionary and plot paths dictionary
    results = {}
    plot_paths = {}
    
    # Create a list to store hyperparameter information
    hyperparameter_info = []
    
    # Generate all parameter combinations
    param_names = list(param_grid.keys())
    param_values = list(param_grid.values())
    param_combinations = list(itertools.product(*param_values))
    
    # Explore all parameter combinations
    for i, combo in enumerate(param_combinations):
        # Create parameter dictionary for this combination
        params = {param_names[j]: combo[j] for j in range(len(param_names))}
        
        # Create key for this parameter combination
        key = "_".join([f"{k}{v}" for k, v in params.items()])
        print(f"\nRunning parameter combination {i+1}/{len(param_combinations)}: {key}")
        
        # Ensure random_state is consistently used
        params_with_random_state = params.copy()
        params_with_random_state['random_state'] = random_state
        
        # Run dimensionality reduction
        result = reduction_method.fit_transform(expr_data, **params_with_random_state)
        
        # Save embedding results
        method_dir = create_parameter_dir(output_dir, reduction_method.get_method_name(), params)
        result_path = os.path.join(method_dir, f"embedding_coordinates.csv")
        result.to_csv(result_path)
        print(f"Saved embedding coordinates to: {result_path}")
        
        # Store hyperparameter information
        param_info = params.copy()
        param_info['embedding_path'] = result_path
        hyperparameter_info.append(param_info)
        
        # Generate and save plots
        paths = plot_and_save_embedding(
            result, metadata,
            method_name=reduction_method.get_method_name(),
            color_by=color_by,
            categorical=categorical,
            params=params,
            output_dir=output_dir
        )
        
        # Store results and paths
        results[key] = {
            'result': result,
            'params': params,
            'embedding_path': result_path
        }
        
        plot_paths[key] = paths
    
    # Create and save hyperparameter dataframe
    hyperparameter_df = pd.DataFrame(hyperparameter_info)
    hyperparameter_df_path = os.path.join(output_dir, f"{method_name}_hyperparameter_map.csv")
    hyperparameter_df.to_csv(hyperparameter_df_path, index=False)
    print(f"Saved hyperparameter mapping to: {hyperparameter_df_path}")
    
    return results, adata_grouped, plot_paths, hyperparameter_df