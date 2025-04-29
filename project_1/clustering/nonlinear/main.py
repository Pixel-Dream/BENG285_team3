import os
import numpy as np
import scanpy as sc
import json
import pandas as pd
from preprocessing import preprocess_metadata
from exploration import run_dimensionality_reduction_with_plots, explore_parameters_with_plots
from visualization import visualize_metadata_distribution

# Example 1: Basic analysis with individual plots
def basic_analysis(adata_path, group_mapping_path, method_name='tsne', color_by=None,
                   params=None, output_dir=None, random_state=42):
    """
    Run a basic dimensionality reduction analysis with plotted results
    
    Parameters:
    -----------
    adata_path : str
        Path to AnnData h5ad file
    group_mapping_path : str
        Path to group mapping JSON file
    method_name : str
        Dimensionality reduction method ('tsne' or 'umap')
    color_by : list of str or None
        List of metadata columns to color by
    params : dict or None
        Method parameters (if None, use defaults)
    output_dir : str or None
        Directory for output plots and files (if None, auto-generate)
    random_state : int or None
        Random state seed for reproducibility
    
    Returns:
    --------
    adata_grouped : AnnData
        AnnData with grouped metadata
    result : DataFrame
        Embedding coordinates
    embedding_path : str
        Path to the saved embedding coordinates CSV
    """
    # Load data
    print(f"Loading AnnData from {adata_path}")
    adata = sc.read_h5ad(adata_path)
    
    # Create output directory if not specified
    if output_dir is None:
        output_dir = f"{method_name.lower()}_plots"
    os.makedirs(output_dir, exist_ok=True)
    
    # Run full analysis with individual plots
    print(f"Running {method_name} analysis with metadata grouping from {group_mapping_path}")
    result, metadata, adata_grouped, plot_paths, embedding_path = run_dimensionality_reduction_with_plots(
        adata, 
        method_name=method_name,
        params=params,
        group_mapping=group_mapping_path,
        output_dir=output_dir,
        random_state=random_state,
        n_top_genes=2000,
        color_by=color_by
    )
    
    # Visualize metadata distributions
    # print("Generating metadata distribution plot")
    # fig = visualize_metadata_distribution(
    #     adata_grouped, 
    #     output_dir=output_dir,
    #     filename="metadata_distributions.png"
    # )
    
    print(f"Analysis complete! Individual plots saved to {output_dir}")
    print(f"Embedding coordinates saved to {embedding_path}")
    
    return adata_grouped, result, embedding_path


# Example 2: Parameter exploration
def parameter_exploration(adata_path, group_mapping_path, method_name='tsne', color_by=None,
                         param_grid=None, output_dir=None, random_state=42):
    """
    Explore dimensionality reduction parameters and generate plots for each configuration
    
    Parameters:
    -----------
    adata_path : str
        Path to AnnData h5ad file
    group_mapping_path : str
        Path to group mapping JSON file
    method_name : str
        Dimensionality reduction method ('tsne' or 'umap')
    color_by : list of str or None
        List of metadata columns to color by
    param_grid : dict or None
        Grid of parameters to explore (if None, use defaults)
    output_dir : str or None
        Directory for output plots and files (if None, auto-generate)
    random_state : int or None
        Random state seed for reproducibility
    
    Returns:
    --------
    results : dict
        Dictionary of results for each parameter combination
    adata_grouped : AnnData
        AnnData with grouped metadata
    hyperparameter_df : DataFrame
        DataFrame with hyperparameters and file paths to saved embeddings
    """
    # Load data
    print(f"Loading AnnData from {adata_path}")
    adata = sc.read_h5ad(adata_path)
    
    # Create output directory if not specified
    if output_dir is None:
        output_dir = f"{method_name.lower()}_parameter_exploration"
    os.makedirs(output_dir, exist_ok=True)
    
    # Run parameter exploration
    print(f"Running {method_name} parameter exploration with metadata grouping from {group_mapping_path}")
    results, adata_grouped, plot_paths, hyperparameter_df = explore_parameters_with_plots(
        adata,
        method_name=method_name,
        param_grid=param_grid,
        group_mapping=group_mapping_path,
        color_by=color_by,
        output_dir=output_dir,
        random_state=random_state
    )
    
    print(f"Parameter exploration complete! Results saved to {output_dir}")
    hyperparameter_path = os.path.join(output_dir, f"{method_name}_hyperparameter_map.csv")
    print(f"Hyperparameter mapping saved to: {hyperparameter_path}")
    
    return results, adata_grouped, hyperparameter_df


# Execute examples if run as script
if __name__ == "__main__":
    # Path to your data files
    ADATA_PATH = "../../data/z-scaled_w_normalized_merged.h5ad"  # Update with your file path
    GROUP_MAPPING_PATH = "../../data/meta_group.json"  # Update with your file path
    
    # Check if data files exist
    if not os.path.exists(ADATA_PATH):
        print(f"Warning: {ADATA_PATH} not found. Please update the path.")
    elif not os.path.exists(GROUP_MAPPING_PATH):
        print(f"Warning: {GROUP_MAPPING_PATH} not found. Please update the path.")
    else:
        # Run basic t-SNE analysis
        print("\n=== Running Basic t-SNE Analysis ===")
        adata_grouped, tsne_result, tsne_path = basic_analysis(
            ADATA_PATH, 
            GROUP_MAPPING_PATH, 
            method_name='tsne',
            output_dir="tsne_plots",
            random_state=123
        )
        
        # Run basic UMAP analysis
        print("\n=== Running Basic UMAP Analysis ===")
        adata_grouped, umap_result, umap_path = basic_analysis(
            ADATA_PATH, 
            GROUP_MAPPING_PATH, 
            method_name='umap',
            output_dir="umap_plots",
            random_state=123
        )
        
       
        
        # t-SNE parameter exploration
        print("\n=== Running t-SNE Parameter Exploration ===")
        tsne_results, adata_grouped, tsne_hyperparameter_df = parameter_exploration(
            ADATA_PATH, 
            GROUP_MAPPING_PATH, 
            method_name='tsne',
            output_dir="tsne_parameter_exploration",
            random_state=123
        )
        
        # UMAP parameter exploration
        print("\n=== Running UMAP Parameter Exploration ===")
        umap_results, adata_grouped, umap_hyperparameter_df = parameter_exploration(
            ADATA_PATH, 
            GROUP_MAPPING_PATH, 
            method_name='umap',
            output_dir="umap_parameter_exploration",
            random_state=123
        )
