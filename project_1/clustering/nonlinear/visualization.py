import os
import numpy as np
import pandas as pd
import matplotlib.pyplot as plt
import seaborn as sns
import re
from pathlib import Path

from utils import sanitize_for_filename, truncate_string, create_parameter_dir

def plot_and_save_embedding(embedding_result, metadata, method_name='Embedding', 
                          color_by=None, categorical=None, params=None,
                          continuous_cmap='viridis', point_size=30, alpha=0.7, 
                          output_dir='plots', max_legend_length=25, max_legend_items=10,
                          figsize=(10, 8)):
    """
    Plot and save embedding results. Always generates an uncolored plot first,
    then generates individual high-quality figures colored by specified metadata.
    
    Parameters:
    -----------
    embedding_result : DataFrame
        DataFrame with embedding coordinates
    metadata : DataFrame
        Metadata with samples as index
    method_name : str
        Name of the dimensionality reduction method
    color_by : str, list, or None
        Metadata column(s) to color points by. If None, only the uncolored plot is generated.
        If a string or list, colored plots are generated in addition to the uncolored one.
    categorical : bool or list
        Whether each color_by column is categorical (used only if color_by is not None).
    params : dict
        Dictionary of parameters used for the method
    continuous_cmap : str
        Colormap for continuous variables
    point_size : int
        Size of scatter plot points
    alpha : float
        Transparency of points
    output_dir : str
        Base output directory
    max_legend_length : int
        Maximum length of each legend entry before truncation
    max_legend_items : int
        Maximum number of legend items to show directly (others grouped)
    figsize : tuple
        Figure size for individual plots
        
    Returns:
    --------
    plot_paths : list
        List of paths to saved plot files (includes the uncolored plot and any colored plots).
    """
    # Set default parameters if not provided
    if params is None:
        params = {}
    
    # Create parameter directory
    param_dir = create_parameter_dir(output_dir, method_name, params)
    
    # Initialize list of plot paths
    plot_paths = []
    
    # Get coordinate column names
    coord_cols = embedding_result.columns
    x_col, y_col = coord_cols[0], coord_cols[1]
    
    # Create short parameter string for filenames
    method_short = method_name.lower().replace('-', '')
    
    # --- Always generate uncolored plot first ---
    fig_uncolored, ax_uncolored = plt.subplots(figsize=figsize, dpi=100)
    
    ax_uncolored.scatter(
        embedding_result[x_col], embedding_result[y_col],
        s=point_size, alpha=alpha, color='gray' # Use a default color
    )
    
    ax_uncolored.set_title(f"{method_name} Visualization (Uncolored)", fontsize=12)
    ax_uncolored.set_xlabel(f'{x_col}', fontsize=10)
    ax_uncolored.set_ylabel(f'{y_col}', fontsize=10)
    
    ax_uncolored.spines['top'].set_visible(False)
    ax_uncolored.spines['right'].set_visible(False)
    ax_uncolored.tick_params(labelsize='small')
    
    uncolored_filename = f"{method_short}_params_uncolored.png"
    uncolored_filepath = os.path.join(param_dir, uncolored_filename)
    
    plt.tight_layout()
    fig_uncolored.savefig(uncolored_filepath, dpi=300, bbox_inches='tight')
    plot_paths.append(uncolored_filepath)
    print(f"Saved uncolored plot to: {uncolored_filepath}")
    plt.close(fig_uncolored)
    # --- End of uncolored plot generation ---

    # --- Generate colored plots if color_by is specified ---
    if color_by is not None and color_by:  # Check if color_by is not None and not empty
        # Handle single metadata column if specified
        if isinstance(color_by, str):
            color_by = [color_by]
            if isinstance(categorical, bool):
                categorical = [categorical]
        elif not isinstance(color_by, list):
             print(f"Warning: Invalid type for color_by: {type(color_by)}. Skipping colored plots.")
             return plot_paths # Return only the uncolored plot path
        
        # Check if categorical list matches color_by list length if provided
        if isinstance(categorical, list) and len(categorical) != len(color_by):
            print(f"Warning: Length of 'categorical' ({len(categorical)}) does not match length of 'color_by' ({len(color_by)}). Using default behavior.")
            categorical = [True] * len(color_by) # Default to categorical or adjust as needed
        elif isinstance(categorical, bool):
            categorical = [categorical] * len(color_by)
        elif not isinstance(categorical, list):
            print(f"Warning: Invalid type for categorical: {type(categorical)}. Assuming all are categorical.")
            categorical = [True] * len(color_by)
            
        # Combine embedding results with metadata for coloring
        data = embedding_result.copy()
        
        # Plot each metadata variable separately
        for i, (col, is_cat) in enumerate(zip(color_by, categorical)):
            # Skip if column not found in metadata
            if col not in metadata.columns:
                print(f"Warning: Column '{col}' not found in metadata")
                continue
            
            # Add metadata column to embedding results
            # Ensure index alignment if metadata index doesn't match embedding_result index
            if not metadata.index.equals(data.index):
                 # Attempt to align based on index
                 try:
                     aligned_metadata_col = metadata.loc[data.index, col]
                 except KeyError:
                     print(f"Warning: Could not align metadata index for column '{col}'. Skipping this column.")
                     continue
                 data[col] = aligned_metadata_col.values
            else:
                data[col] = metadata[col].values
            
            # Create a clean, short column name for the plot title and filename
            clean_col = col.split('.')[-1] if '.' in col else col
            col_title = clean_col.replace('_', ' ').title()
            
            # Create figure for colored plot
            fig_colored, ax_colored = plt.subplots(figsize=figsize, dpi=100)
            
            if is_cat:  # Categorical variables
                # Convert to string to handle numeric categories safely
                data[col] = data[col].astype(str)
                
                # Handle potential NaN values after conversion (though less likely)
                data[col] = data[col].fillna('NaN')
                
                unique_categories = data[col].unique()
                n_categories = len(unique_categories)
                
                if n_categories > max_legend_items:
                    top_categories = data[col].value_counts().nlargest(max_legend_items-1).index.tolist()
                    def map_category(x): return x if x in top_categories else 'Other'
                    data['mapped_category'] = data[col].apply(map_category)
                    other_count = (data['mapped_category'] == 'Other').sum()
                    plot_col = 'mapped_category'
                    print(f"Column '{col}' has {n_categories} categories; showing top {max_legend_items-1} plus 'Other' ({other_count} items)")
                else:
                    plot_col = col
                
                categories = sorted(data[plot_col].unique())
                palette = sns.color_palette('husl', n_colors=len(categories))
                
                for j, category in enumerate(categories):
                    truncated_category = truncate_string(category, max_legend_length)
                    mask = data[plot_col] == category
                    ax_colored.scatter(
                        data.loc[mask, x_col], data.loc[mask, y_col],
                        s=point_size, alpha=alpha, label=truncated_category, color=palette[j]
                    )
                
                legend = ax_colored.legend(
                    title=col_title, loc='best', fontsize='small',
                    markerscale=0.8, framealpha=0.9
                )
                legend.get_title().set_fontsize('small')
                
            else:  # Continuous variables
                # Attempt conversion to numeric, handle NaNs
                numeric_col = pd.to_numeric(data[col], errors='coerce')
                
                if numeric_col.isna().all():
                    ax_colored.text(0.5, 0.5, f"All values in '{col}' are non-numeric or NaN", 
                            ha='center', va='center', transform=ax_colored.transAxes)
                else:
                    # Plot continuous data, NaNs will be ignored by scatter
                    scatter = ax_colored.scatter(
                        data[x_col], data[y_col],
                        s=point_size, alpha=alpha, c=numeric_col, cmap=continuous_cmap
                    )
                    cbar = plt.colorbar(scatter, ax=ax_colored, label=col_title)
                    cbar.ax.tick_params(labelsize='small')
            
            # Set title and labels for colored plot
            ax_colored.set_title(f"{method_name} - {col_title}", fontsize=12)
            ax_colored.set_xlabel(f'{x_col}', fontsize=10)
            ax_colored.set_ylabel(f'{y_col}', fontsize=10)
            
            # Clean up plot appearance
            ax_colored.spines['top'].set_visible(False)
            ax_colored.spines['right'].set_visible(False)
            ax_colored.tick_params(labelsize='small')
            
            # Create clean filename for colored plot
            clean_filename = f"{method_short}_params_{sanitize_for_filename(clean_col)}.png"
            filepath = os.path.join(param_dir, clean_filename)
            
            # Save figure with high DPI
            plt.tight_layout()
            fig_colored.savefig(filepath, dpi=300, bbox_inches='tight')
            plot_paths.append(filepath)
            
            print(f"Saved colored plot to: {filepath}")
            plt.close(fig_colored) # Close the figure for this colored plot
            
            # Remove temporary columns if created
            if 'mapped_category' in data.columns:
                data = data.drop(columns=['mapped_category'])
            if col not in embedding_result.columns: # Avoid removing original embedding columns if name overlaps
                 data = data.drop(columns=[col])
                 
    # --- End of colored plots generation ---
    
    return plot_paths


def visualize_metadata_distribution(adata, grouped_columns=None, figsize=(15, 10), 
                                   output_dir=None, filename="metadata_distribution.png"):
    """
    Visualize the distribution of metadata categories
    
    Parameters:
    -----------
    adata : AnnData
        AnnData object with metadata
    grouped_columns : list or None
        List of grouped columns to visualize; if None, visualize all columns ending with '_grouped'
    figsize : tuple
        Figure size
    output_dir : str or None
        Output directory; if None, figure is not saved
    filename : str
        Filename for saved figure
        
    Returns:
    --------
    fig : Figure
        Matplotlib figure
    """
    import matplotlib.pyplot as plt
    import seaborn as sns
    
    # If no columns specified, find all grouped columns
    if grouped_columns is None:
        grouped_columns = [col for col in adata.obs.columns if col.endswith('_grouped')]
    
    # Determine number of plots needed
    n_cols = 2
    n_rows = (len(grouped_columns) + n_cols - 1) // n_cols
    
    # Create figure and axes
    fig, axes = plt.subplots(n_rows, n_cols, figsize=figsize)
    
    # Flatten axes for easy iteration
    if n_rows * n_cols > 1:
        axes = axes.flatten()
    else:
        axes = [axes]
    
    # Plot each metadata distribution
    for i, col in enumerate(grouped_columns):
        if i < len(axes):
            # Get value counts
            value_counts = adata.obs[col].value_counts()
            
            # Create bar chart
            sns.barplot(x=value_counts.index, y=value_counts.values, ax=axes[i])
            
            # Format axis
            axes[i].set_title(f'Distribution of {col}')
            axes[i].set_xlabel('')
            axes[i].set_ylabel('Count')
            
            # Rotate x-axis labels if many categories
            if len(value_counts) > 4:
                axes[i].set_xticklabels(axes[i].get_xticklabels(), rotation=45, ha='right')
                plt.setp(axes[i].xaxis.get_majorticklabels(), ha='right')
    
    # Hide unused subplots
    for i in range(len(grouped_columns), len(axes)):
        axes[i].set_visible(False)
    
    plt.tight_layout()
    
    # Save figure if output_dir provided
    if output_dir is not None:
        os.makedirs(output_dir, exist_ok=True)
        filepath = os.path.join(output_dir, filename)
        fig.savefig(filepath, dpi=300, bbox_inches='tight')
        print(f"Saved metadata distribution to: {filepath}")
    
    return fig