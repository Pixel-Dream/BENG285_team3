#%%
import os
import matplotlib.pyplot as plt
import matplotlib.cm as cm
from matplotlib.colors import to_rgba
import pandas as pd
import numpy as np
import scanpy as sc

import sys
sys.path.append("/Users/xbh0403/Desktop/25SP/BENG285/reanalyze/clustering")
from preprocessing import preprocess_anndata, preprocess_metadata

adata_haowen = sc.read_h5ad("/Users/xbh0403/Desktop/25SP/BENG285/tsne_analysis/z_scaled_w_normalized_merged.h5ad")

adata_haowen = preprocess_metadata(adata_haowen, "/Users/xbh0403/Desktop/25SP/BENG285/tsne_analysis/meta_group.json", inplace=False)
expr_data, metadata = preprocess_anndata(adata_haowen)

groups = [
    "demographic.country_of_residence_at_enrollment_grouped",
    "demographic.gender_grouped",
    "demographic.race_grouped",
    "demographic.vital_status_grouped",
    "diagnoses.ajcc_pathologic_m_grouped",
    "diagnoses.ajcc_pathologic_n_grouped",
    "diagnoses.ajcc_pathologic_stage_grouped",
    "diagnoses.ajcc_pathologic_t_grouped",
    "diagnoses.classification_of_tumor_grouped",
    "diagnoses.diagnosis_is_primary_disease_grouped",
    "diagnoses.laterality_grouped",
    "diagnoses.primary_diagnosis_grouped",
    "diagnoses.prior_malignancy_grouped",
    "diagnoses.prior_treatment_grouped",
    "diagnoses.sites_of_involvement_grouped",
    "diagnoses.tissue_or_organ_of_origin_grouped",
    "follow_ups.disease_response_grouped",
    "cases.disease_type_grouped",
]

#%%

def plot_pie_chart(df, column_name, title=None, figsize=(10, 8),
                   colors=None, startangle=90,
                   threshold=0.03, explode=None, show_legend=True,
                   style='publication', font_family='sans-serif', 
                   shadow=False, show_values=True, text_color='auto',
                   edge_width=0.8, edge_color='white', label_distance=1.15,
                   save_path_dir="/Users/xbh0403/Desktop/25SP/BENG285/figures/others", dpi=300, save_format='png', transparent=False):
    """
    Creates a publication-quality pie chart from a column in a dataframe.
    
    Parameters:
    -----------
    df : pandas.DataFrame
        The dataframe containing the data
    column_name : str
        The name of the column to use for the pie chart
    title : str, optional
        The title of the pie chart. If None, uses "Distribution of {column_name}"
    figsize : tuple, optional
        The size of the figure (width, height)
    colors : list or str, optional
        List of colors for the pie slices or name of a colormap. Options include:
        'publication', 'pastel', 'bright', 'muted', 'colorblind' or any matplotlib colormap
    startangle : int, optional
        The angle to start the pie chart from
    threshold : float, optional
        Group categories with percentage less than threshold into "Other"
    explode : list or float, optional
        List of values to "explode" each slice by or a single value for all slices
    show_legend : bool, optional
        Whether to show the legend
    style : str, optional
        Overall style of the chart. Options: 'publication', 'modern', 'minimal', 'classic'
    font_family : str, optional
        Font family for all text elements. Common choices: 'sans-serif', 'serif'
    shadow : bool, optional
        Whether to add shadow effect (not recommended for publication)
    show_values : bool, optional
        Whether to show values and percentages on the pie
    text_color : str, optional
        Color for text on pie slices. 'auto' will choose black or white based on background
    edge_width : float, optional
        Width of the edge between slices
    edge_color : str, optional
        Color of the edge between slices
    label_distance : float, optional
        Distance of labels from the center (1.0 is at the edge)
    save_path : str, optional
        Path to save the figure. If None, figure is not saved
    dpi : int, optional
        Resolution of the saved figure in dots per inch
    save_format : str, optional
        Format to save the figure in ('png', 'pdf', 'svg', 'jpg', etc.)
    transparent : bool, optional
        Whether to save the figure with a transparent background
        
    Returns:
    --------
    fig, ax : tuple
        The figure and axis objects for further customization if needed
    """
    save_path = os.path.join(save_path_dir, f"{column_name}.{save_format}")
    
    # Validate inputs
    if column_name not in df.columns:
        raise ValueError(f"Column '{column_name}' not found in dataframe. Available columns: {df.columns.tolist()}")

    # Apply the chosen style
    with plt.style.context(['seaborn-v0_8-whitegrid', 'seaborn-v0_8-paper']):
        
        # Set specific style parameters based on style choice
        if style == 'publication':
            plt.rcParams.update({
                'font.family': font_family,
                'font.size': 11,
                'axes.labelsize': 12,
                'axes.titlesize': 14,
                'xtick.labelsize': 10,
                'ytick.labelsize': 10,
                'legend.fontsize': 10,
                'figure.titlesize': 16
            })
        elif style == 'modern':
            plt.rcParams.update({
                'font.family': font_family,
                'font.size': 12,
                'axes.labelsize': 13,
                'axes.titlesize': 15,
                'xtick.labelsize': 11,
                'ytick.labelsize': 11,
                'legend.fontsize': 11,
                'figure.titlesize': 17
            })
        elif style == 'minimal':
            plt.rcParams.update({
                'font.family': font_family,
                'font.size': 10,
                'axes.labelsize': 11,
                'axes.titlesize': 13,
                'xtick.labelsize': 9,
                'ytick.labelsize': 9,
                'legend.fontsize': 9,
                'figure.titlesize': 14
            })
        elif style == 'classic':
            plt.rcParams.update({
                'font.family': 'serif',
                'font.size': 11,
                'axes.labelsize': 12,
                'axes.titlesize': 14,
                'xtick.labelsize': 10,
                'ytick.labelsize': 10,
                'legend.fontsize': 10,
                'figure.titlesize': 16
            })

        # Get value counts
        value_counts = df[column_name].value_counts()
        total = value_counts.sum()

        # Calculate percentages for threshold check
        percentages = value_counts / total

        # Group small categories into "Other" if they're below the threshold
        if threshold > 0:
            small_categories = percentages[percentages < threshold].index
            # Check if there are small categories to group
            if len(small_categories) > 1:  # Need at least 2 to group
                # Create a copy to avoid modification warnings
                value_counts_modified = value_counts.copy()

                # Sum the small categories and create an "Other" category
                others_sum = value_counts_modified[small_categories].sum()
                for cat in small_categories:
                    value_counts_modified = value_counts_modified.drop(cat)

                if others_sum > 0:  # Only add "Other" if there's a value
                    value_counts_modified["Other"] = others_sum

                value_counts = value_counts_modified  # Update value_counts for plotting

        # Create plot
        fig, ax = plt.subplots(figsize=figsize, facecolor='white')
        
        # Set the figure background to white
        fig.patch.set_facecolor('white')
        ax.set_facecolor('white')

        # Set up colors based on provided parameter or style
        num_slices = len(value_counts)
        
        # Define color schemes suited for publications
        color_schemes = {
            'publication': ['#4878D0', '#EE854A', '#6ACC64', '#D65F5F', '#956CB4', 
                          '#8C613C', '#DC7EC0', '#82C6E2', '#D5BB67', '#82C6E2'],
            'pastel': ['#c6dbef', '#fdd0a2', '#c7e9c0', '#fcbba1', '#d4b9da', 
                      '#e6baad', '#f2b6d2', '#b3e0dc', '#f3d9a2', '#c4c1e0'],
            'bright': ['#1f77b4', '#ff7f0e', '#2ca02c', '#d62728', '#9467bd', 
                      '#8c564b', '#e377c2', '#7f7f7f', '#bcbd22', '#17becf'],
            'muted': ['#4878D0', '#D65F5F', '#6ACC64', '#956CB4', '#8C613C', 
                     '#DC7EC0', '#82C6E2', '#D5BB67', '#C4C1E0', '#BBB0AC'],
            'colorblind': ['#0173B2', '#DE8F05', '#029E73', '#D55E00', '#CC78BC', 
                          '#CA9161', '#FBAFE4', '#949494', '#ECE133', '#56B4E9']
        }

        if colors is None:
            # Use publication scheme as default
            colors = color_schemes.get(style, color_schemes['publication'])
            # Repeat colors if needed
            colors = [colors[i % len(colors)] for i in range(num_slices)]
        elif isinstance(colors, str):
            # If a string is provided, check if it's one of our schemes
            if colors in color_schemes:
                colors = [color_schemes[colors][i % len(color_schemes[colors])] for i in range(num_slices)]
            else:
                # Otherwise, treat it as a matplotlib colormap name
                try:
                    cmap = cm.get_cmap(colors)
                    colors = [cmap(i/num_slices) for i in range(num_slices)]
                except:
                    print(f"Warning: Unknown colormap '{colors}'. Using default.")
                    colors = [color_schemes['publication'][i % len(color_schemes['publication'])] for i in range(num_slices)]

        # Set up explode parameter
        if explode is None:
            # No explode by default for publication style
            explode = [0] * num_slices
        elif isinstance(explode, (int, float)):
            # If a single value is provided, use it for all slices
            explode = [explode] * num_slices
        elif len(explode) != num_slices:
            print(f"Warning: Explode length mismatch. Using default explode. Expected {num_slices}, got {len(explode)}")
            explode = [0] * num_slices

        # Define wedge properties for better separation
        wedgeprops = {
            'linewidth': edge_width, 
            'edgecolor': edge_color,
            'antialiased': True
        }

        # Function to determine text color based on background
        def get_text_color(bg_color):
            if text_color != 'auto':
                return text_color
                
            # Convert color to RGB
            if isinstance(bg_color, str):
                rgba = to_rgba(bg_color)
            else:
                rgba = bg_color
                
            # Calculate luminance
            r, g, b = rgba[:3]
            luminance = 0.299 * r + 0.587 * g + 0.114 * b
            
            # Return black for light backgrounds, white for dark
            return 'black' if luminance > 0.5 else 'white'

        # Function to display both count and percentage
        def make_autopct(values):
            def my_autopct(pct):
                if not show_values:
                    return ''
                    
                total_val = sum(values)
                val = int(round(pct * total_val / 100.0))
                
                # Don't display label for very small slices to avoid clutter
                if pct < 2:
                    return ''
                elif pct < 5:
                    return f'{pct:.1f}%'
                else:
                    return f'{val:,}\n({pct:.1f}%)'
            return my_autopct

        # Use the custom autopct function
        autopct_func = make_autopct(value_counts.values)

        # Create pie chart
        wedges, texts, autotexts = ax.pie(
            value_counts,
            labels=None,  # No labels on the pie directly
            autopct=autopct_func,
            startangle=startangle,
            explode=explode,
            colors=colors,
            shadow=shadow,
            wedgeprops=wedgeprops,
            pctdistance=0.65,  # Move percentages closer to center
            textprops={'fontsize': 10, 'fontweight': 'medium'},
            radius=0.9  # Slightly smaller pie to leave room for labels/legend
        )

        # Adjust text colors for readability
        for i, autotext in enumerate(autotexts):
            # Get the corresponding pie slice color
            bg_color = colors[i % len(colors)]
            # Set appropriate text color
            autotext.set_color(get_text_color(bg_color))
            autotext.set_fontsize(10)
            autotext.set_fontweight('medium')
            autotext.set_va('center')

        # Equal aspect ratio ensures the pie chart is circular
        ax.axis('equal')

        # Add title
        if title is None:
            title = f"Distribution of {column_name}"
            
        # Adjust title style based on chosen overall style
        if style == 'publication':
            plt.title(title, fontsize=14, pad=20, fontweight='bold')
        elif style == 'modern':
            plt.title(title, fontsize=15, pad=20, fontweight='semibold')
        elif style == 'minimal':
            plt.title(title, fontsize=13, pad=15, fontweight='medium')
        elif style == 'classic':
            plt.title(title, fontsize=14, pad=20, fontweight='medium', 
                     fontfamily='serif')

        # Add legend if requested
        if show_legend:
            # Format legend with sorted values (largest first)
            sorted_items = value_counts.sort_values(ascending=False)

            # --- Choose Legend Label Format ---
            # Option 1: Label (Count, Percentage) - tends to be informative
            legend_labels = [f'{label} ({sorted_items[label]:,}, {percentages[label]:.1%})'
                             for label in sorted_items.index]
            # Option 2: Label (Percentage) - cleaner if counts are not essential
            # legend_labels = [f'{label} ({percentages[label]:.1%})'
            #                  for label in sorted_items.index]
            # Option 3: Label (Count)
            # legend_labels = [f'{label} ({sorted_items[label]:,})'
            #                  for label in sorted_items.index]
            # --- End Label Format Choice ---

            # Get corresponding wedges in the same order
            sorted_indices = [list(value_counts.index).index(idx) for idx in sorted_items.index]
            sorted_wedges = [wedges[i] for i in sorted_indices]

            # --- New Legend Placement Logic ---
            legend = plt.legend(
                sorted_wedges,
                legend_labels,
                title=column_name.replace('_', ' ').title(), # Use cleaned column name as title
                loc='lower center', # Position below the center of the plot area
                bbox_to_anchor=(0.5, -0.15), # Adjust vertical position below the axes (tweak if needed)
                frameon=False, # No frame for cleaner look
                ncol=3, # Max 3 items per row, wraps automatically
                fontsize=plt.rcParams.get('legend.fontsize', 10), # Use style fontsize, default 10
                title_fontsize=plt.rcParams.get('legend.fontsize', 10) # Title same size as items
            )
            legend.get_title().set_fontweight('bold')
            # --- End New Legend Placement Logic ---

        # Adjust layout to prevent labels/legend overlapping
        # Use tight_layout, rect reserves space: [left, bottom, right, top]
        try:
            # Add padding at the bottom (e.g., 0.1 means 10% of fig height)
            # Adjust the second value (bottom) if legend is clipped or too far
            fig.tight_layout(rect=[0, 0.1, 1, 0.95]) # Leave 10% at bottom, 5% at top for title
        except ValueError:
             print("Warning: tight_layout failed. Legend might be clipped or overlap title.")
             # As an alternative, try subplots_adjust if tight_layout consistently fails
             # fig.subplots_adjust(bottom=0.2, top=0.9)

        # For minimal style, remove spines
        if style == 'minimal':
            ax.spines['top'].set_visible(False)
            ax.spines['right'].set_visible(False)
            ax.spines['bottom'].set_visible(False)
            ax.spines['left'].set_visible(False)
            
        # Save the figure if a save path is provided
        if save_path is not None:
            # If save_format is not in the save_path, append it
            if not save_path.lower().endswith(f'.{save_format.lower()}'):
                save_path = f"{save_path}.{save_format.lower()}"
                
            # Create directory if it doesn't exist
            save_dir = os.path.dirname(save_path)
            if save_dir and not os.path.exists(save_dir):
                os.makedirs(save_dir)
                
            # Save the figure
            plt.savefig(
                save_path,
                dpi=dpi,
                bbox_inches='tight',
                transparent=transparent,
                format=save_format
            )
            print(f"Figure saved to: {save_path}")

        return fig, ax

#%%

for group in groups:
    fig, ax = plot_pie_chart(metadata, group, threshold=0)
    plt.show()


#%%

import os
import numpy as np
import matplotlib.pyplot as plt
from PIL import Image
from matplotlib.gridspec import GridSpec

def create_image_grid(image_folder, output_path=None, grid_size=None, 
                      figsize=None, titles=None, spacing=0.1, 
                      dpi=300, border_color='white', border_size=5,
                      bg_color='white', method='matplotlib'):
    """
    Creates a grid of images from all image files in a folder and saves the combined image.
    
    Parameters:
    -----------
    image_folder : str
        Path to the folder containing image files
    output_path : str, optional
        Path to save the combined image. If None, uses 'image_grid.png'
    grid_size : tuple, optional
        Size of the grid (rows, cols). If None, calculated automatically
    figsize : tuple, optional
        Size of the output figure in inches. If None, calculated based on grid size
    titles : list or bool, optional
        List of titles for each image, or True to use filenames, or False for no titles
    spacing : float, optional
        Spacing between images as a fraction of figure size
    dpi : int, optional
        Resolution for saved image
    border_color : str, optional
        Color of the border around each image
    border_size : int, optional
        Size of the border in pixels
    bg_color : str, optional
        Background color for the grid
    method : str, optional
        Method to use for creating the grid ('matplotlib' or 'pil')
        
    Returns:
    --------
    fig : matplotlib.figure.Figure or PIL.Image
        The figure object with the grid of images
    """
    # Get all image files in the folder
    valid_extensions = ['.jpg', '.jpeg', '.png', '.bmp', '.gif', '.tiff', '.tif']
    image_files = [f for f in os.listdir(image_folder) if 
                  os.path.isfile(os.path.join(image_folder, f)) and
                  os.path.splitext(f.lower())[1] in valid_extensions]
    
    if not image_files:
        raise ValueError(f"No image files found in {image_folder}")
    
    # Sort image files for consistent ordering
    image_files.sort()
    
    # Set default output path if not provided
    if output_path is None:
        output_path = 'image_grid.png'
    
    # Ensure the file extension is there
    if '.' not in output_path:
        output_path = f"{output_path}.png"
    
    # Choose method
    if method.lower() == 'pil':
        return create_image_montage_pil(
            image_folder, image_files, output_path, grid_size, 
            padding=int(spacing * 100), bg_color=bg_color, border_size=border_size, 
            border_color=border_color
        )
    else:  # Default to matplotlib
        return create_image_grid_matplotlib(
            image_folder, image_files, output_path, grid_size, figsize,
            titles, spacing, dpi, border_color, border_size, bg_color
        )


def create_image_grid_matplotlib(image_folder, image_files, output_path, 
                                grid_size=None, figsize=None, titles=None, 
                                spacing=0.1, dpi=300, border_color='white', 
                                border_size=5, bg_color='white'):
    """
    Creates a grid of images using matplotlib (better for annotations and titles).
    """
    # Determine grid size if not specified
    n_images = len(image_files)
    if grid_size is None:
        # Calculate a reasonable grid size
        cols = int(np.ceil(np.sqrt(n_images)))
        rows = int(np.ceil(n_images / cols))
        grid_size = (rows, cols)
    else:
        rows, cols = grid_size
        # Ensure the grid is large enough
        if rows * cols < n_images:
            print(f"Warning: Grid size {grid_size} too small for {n_images} images. Adjusting...")
            cols = int(np.ceil(n_images / rows))
            grid_size = (rows, cols)
    
    # Determine figure size if not specified
    if figsize is None:
        # Default to 3 inches per image, adjusted by the number of columns
        figsize = (3 * cols, 3 * rows)
    
    # Create figure and gridspec
    fig = plt.figure(figsize=figsize, dpi=dpi, facecolor=bg_color)
    gs = GridSpec(rows, cols, figure=fig, wspace=spacing, hspace=spacing)
    
    # Process titles
    if titles is True:
        titles = [os.path.splitext(os.path.basename(f))[0] for f in image_files]
    elif titles is None or titles is False:
        titles = [None] * n_images
    
    # Add each image to the grid
    for i, img_file in enumerate(image_files):
        if i >= rows * cols:
            print(f"Warning: Only showing {rows * cols} out of {n_images} images")
            break
        
        # Calculate row and column
        row, col = i // cols, i % cols
        
        # Create subplot
        ax = fig.add_subplot(gs[row, col])
        
        # Load image with PIL
        img_path = os.path.join(image_folder, img_file)
        try:
            img = Image.open(img_path)
            
            # Add border if specified
            if border_size > 0:
                # Create slightly larger image with border color
                if img.mode != 'RGB':
                    img = img.convert('RGB')
                    
                border_img = Image.new('RGB', 
                                      (img.width + 2*border_size, img.height + 2*border_size), 
                                      border_color)
                # Paste original image in center
                border_img.paste(img, (border_size, border_size))
                img = border_img
            
            # Display image
            ax.imshow(np.array(img))
            
            # Add title if available
            if i < len(titles) and titles[i]:
                ax.set_title(titles[i], fontsize=10)
            
            # Remove axes
            ax.axis('off')
            
        except Exception as e:
            print(f"Error loading image {img_file}: {e}")
            ax.text(0.5, 0.5, f"Error loading\n{img_file}", 
                   ha='center', va='center', transform=ax.transAxes)
            ax.axis('off')
    
    # Hide any unused subplots
    for i in range(n_images, rows * cols):
        row, col = i // cols, i % cols
        ax = fig.add_subplot(gs[row, col])
        ax.axis('off')
    
    # Adjust layout
    plt.tight_layout()
    
    # Save figure
    plt.savefig(output_path, dpi=dpi, bbox_inches='tight', facecolor=bg_color)
    print(f"Image grid saved to {output_path}")
    
    return fig


def create_image_montage_pil(image_folder, image_files, output_path, 
                          grid_size=None, image_size=None, padding=10, 
                          bg_color='white', border_size=0, border_color='white'):
    """
    Creates a montage of images using PIL (better for large numbers of images and precise layout).
    """
    # Load all images
    images = []
    for img_file in image_files:
        try:
            img = Image.open(os.path.join(image_folder, img_file))
            if img.mode != 'RGB':
                img = img.convert('RGB')
            images.append(img)
        except Exception as e:
            print(f"Error loading image {img_file}: {e}")
    
    if not images:
        raise ValueError("No valid images could be loaded")
    
    # Determine grid size if not specified
    n_images = len(images)
    if grid_size is None:
        # Calculate a reasonable grid size
        cols = int(np.ceil(np.sqrt(n_images)))
        rows = int(np.ceil(n_images / cols))
        grid_size = (rows, cols)
    else:
        rows, cols = grid_size
        # Ensure the grid is large enough
        if rows * cols < n_images:
            print(f"Warning: Grid size {grid_size} too small for {n_images} images. Adjusting...")
            cols = int(np.ceil(n_images / rows))
            grid_size = (rows, cols)
    
    # Determine uniform image size if not specified
    if image_size is None:
        # Find the maximum dimensions
        max_width = max(img.width for img in images)
        max_height = max(img.height for img in images)
        image_size = (max_width, max_height)
    
    # Add borders if specified
    if border_size > 0:
        for i in range(len(images)):
            img = images[i]
            # Create slightly larger image with border color
            border_img = Image.new('RGB', 
                                  (img.width + 2*border_size, img.height + 2*border_size), 
                                  border_color)
            # Paste original image in center
            border_img.paste(img, (border_size, border_size))
            images[i] = border_img
            
        # Update image_size to include borders
        image_size = (image_size[0] + 2*border_size, image_size[1] + 2*border_size)
    
    # Calculate the total size of the montage
    montage_width = cols * image_size[0] + (cols - 1) * padding
    montage_height = rows * image_size[1] + (rows - 1) * padding
    
    # Create a new image for the montage
    montage = Image.new('RGB', (montage_width, montage_height), bg_color)
    
    # Place images in the montage
    for i, img in enumerate(images):
        if i >= rows * cols:
            print(f"Warning: Only showing {rows * cols} out of {n_images} images")
            break
        
        # Calculate position
        row, col = i // cols, i % cols
        x = col * (image_size[0] + padding)
        y = row * (image_size[1] + padding)
        
        # Resize if needed
        if img.size != image_size:
            # Create a new blank image of the desired size
            resized = Image.new('RGB', image_size, bg_color)
            
            # Calculate position to center the original image
            paste_x = (image_size[0] - img.width) // 2
            paste_y = (image_size[1] - img.height) // 2
            
            # Paste the original image onto the blank image
            resized.paste(img, (paste_x, paste_y))
            img = resized
        
        # Place image in montage
        montage.paste(img, (x, y))
    
    # Save the montage
    montage.save(output_path, quality=95)
    print(f"Image montage saved to {output_path}")
    
    return montage


# Function to create a grid with automatic thumbnail generation
def create_thumbnail_grid(image_folder, output_path=None, grid_size=None,
                        thumb_size=(300, 300), spacing=10, bg_color='white',
                        border_size=5, border_color='white', titles=False):
    """
    Creates a grid of image thumbnails from a folder.
    
    Parameters:
    -----------
    image_folder : str
        Path to the folder containing image files
    output_path : str, optional
        Path to save the combined image
    grid_size : tuple, optional
        Size of the grid (rows, cols)
    thumb_size : tuple, optional
        Size of each thumbnail (width, height)
    spacing : int, optional
        Spacing between thumbnails in pixels
    bg_color : str, optional
        Background color
    border_size : int, optional
        Size of the border around each thumbnail
    border_color : str, optional
        Color of the border
    titles : bool, optional
        Whether to add filenames as titles
        
    Returns:
    --------
    PIL.Image
        The combined image
    """
    # Get all image files
    valid_extensions = ['.jpg', '.jpeg', '.png', '.bmp', '.gif', '.tiff', '.tif']
    image_files = [f for f in os.listdir(image_folder) if 
                  os.path.isfile(os.path.join(image_folder, f)) and
                  os.path.splitext(f.lower())[1] in valid_extensions]
    
    if not image_files:
        raise ValueError(f"No image files found in {image_folder}")
    
    # Sort image files
    image_files.sort()
    
    # Create thumbnails
    thumbnails = []
    for img_file in image_files:
        try:
            # Open image
            img_path = os.path.join(image_folder, img_file)
            img = Image.open(img_path)
            
            # Convert to RGB if needed
            if img.mode != 'RGB':
                img = img.convert('RGB')
            
            # Calculate resize dimensions while preserving aspect ratio
            img_width, img_height = img.size
            ratio = min(thumb_size[0] / img_width, thumb_size[1] / img_height)
            new_size = (int(img_width * ratio), int(img_height * ratio))
            
            # Resize image
            img = img.resize(new_size, Image.LANCZOS)
            
            # Create a blank thumbnail of the exact thumb_size
            thumb = Image.new('RGB', thumb_size, bg_color)
            
            # Calculate position to center the image in the thumbnail
            x = (thumb_size[0] - new_size[0]) // 2
            y = (thumb_size[1] - new_size[1]) // 2
            
            # Paste the resized image onto the blank thumbnail
            thumb.paste(img, (x, y))
            
            # Add border if specified
            if border_size > 0:
                # Create a larger image with the border
                bordered = Image.new('RGB', 
                                    (thumb_size[0] + 2*border_size, 
                                     thumb_size[1] + 2*border_size), 
                                    border_color)
                bordered.paste(thumb, (border_size, border_size))
                thumb = bordered
            
            thumbnails.append((thumb, img_file))
            
        except Exception as e:
            print(f"Error creating thumbnail for {img_file}: {e}")
    
    # Determine grid dimensions if not specified
    n_images = len(thumbnails)
    if grid_size is None:
        cols = int(np.ceil(np.sqrt(n_images)))
        rows = int(np.ceil(n_images / cols))
    else:
        rows, cols = grid_size
        # Ensure the grid is large enough
        if rows * cols < n_images:
            cols = int(np.ceil(n_images / rows))
    
    # Calculate the thumbnail dimensions including border
    if border_size > 0:
        thumb_width = thumb_size[0] + 2*border_size
        thumb_height = thumb_size[1] + 2*border_size
    else:
        thumb_width, thumb_height = thumb_size
    
    # Add space for titles if needed
    title_height = 30 if titles else 0
    
    # Calculate total grid dimensions
    grid_width = cols * thumb_width + (cols - 1) * spacing
    grid_height = rows * (thumb_height + title_height) + (rows - 1) * spacing
    
    # Create grid image
    grid_img = Image.new('RGB', (grid_width, grid_height), bg_color)
    
    # Place thumbnails in the grid
    for i, (thumb, filename) in enumerate(thumbnails):
        if i >= rows * cols:
            break
        
        row, col = i // cols, i % cols
        x = col * (thumb_width + spacing)
        y = row * (thumb_height + title_height + spacing)
        
        # Paste thumbnail
        grid_img.paste(thumb, (x, y))
        
        # Add filename as title if requested
        if titles:
            # Create a blank image for the title
            title_img = Image.new('RGB', (thumb_width, title_height), bg_color)
            # Create a drawing context
            from PIL import ImageDraw, ImageFont
            draw = ImageDraw.Draw(title_img)
            # Try to use a standard font, or default to the PIL default
            try:
                font = ImageFont.truetype("arial.ttf", 12)
            except IOError:
                font = ImageFont.load_default()
            
            # Truncate filename if too long
            title_text = os.path.splitext(filename)[0]
            if len(title_text) > 25:
                title_text = title_text[:22] + "..."
            
            # Calculate text position to center it
            text_width = draw.textlength(title_text, font=font)
            text_x = (thumb_width - text_width) // 2
            
            # Draw the text
            draw.text((text_x, 5), title_text, fill='black', font=font)
            
            # Paste title below the thumbnail
            grid_img.paste(title_img, (x, y + thumb_height))
    
    # Set default output path if not provided
    if output_path is None:
        output_path = 'thumbnail_grid.png'
    
    # Ensure the file extension is there
    if '.' not in output_path:
        output_path = f"{output_path}.png"
    
    # Save the grid
    grid_img.save(output_path, quality=95)
    print(f"Thumbnail grid saved to {output_path}")
    
    return grid_img


# Example folder path (replace with your actual folder path)
folder_path = "/Users/xbh0403/Desktop/25SP/BENG285/figures/others"

# Create a grid with matplotlib (good for adding titles and annotations)
fig = create_image_grid(
    folder_path,
    output_path="grid_matplotlib_others.png",
    grid_size=(3, 6),  # Auto-calculate
    titles=False,     # Use filenames as titles
    spacing=0.1,
    border_size=5,
    border_color='white',
    method='matplotlib'
)
