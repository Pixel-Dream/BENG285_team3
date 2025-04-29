import os
import re
import pandas as pd

def sanitize_for_filename(s):
    """
    Convert a string to a valid filename by replacing invalid characters
    
    Parameters:
    -----------
    s : str
        String to sanitize
        
    Returns:
    --------
    sanitized : str
        Sanitized string suitable for a filename
    """
    # Replace any character that's not alphanumeric, underscore, or dash with underscore
    return re.sub(r'[^\w\-]', '_', str(s))


def truncate_string(s, max_length=25):
    """
    Truncate a string to a maximum length and add ellipsis if necessary
    
    Parameters:
    -----------
    s : str
        String to truncate
    max_length : int
        Maximum length of the truncated string
        
    Returns:
    --------
    truncated : str
        Truncated string
    """
    if isinstance(s, str) and len(s) > max_length:
        return s[:max_length-3] + '...'
    return s


def create_parameter_dir(base_dir, method_name, params):
    """
    Create and return the path to a parameter-specific directory
    
    Parameters:
    -----------
    base_dir : str
        Base output directory
    method_name : str
        Name of the dimensionality reduction method (e.g., 't-SNE', 'UMAP')
    params : dict
        Dictionary of parameters
        
    Returns:
    --------
    param_dir : str
        Path to parameter-specific directory
    """
    # Create parameter string for directory naming
    if method_name.lower() == 't-sne':
        # t-SNE parameters
        p = params.get('perplexity', 30)
        lr = params.get('learning_rate', 200)
        it = params.get('n_iter', 1000)
        ee = params.get('early_exaggeration', 12)
        init = params.get('init', 'pca')
        metric = params.get('metric', 'euclidean')
        param_str = f"tsne_p{p}_lr{lr}_it{it}_ee{ee}_init{init}_met{metric}"
    elif method_name.lower() == 'umap':
        # UMAP parameters
        n = params.get('n_neighbors', 15)
        md = params.get('min_dist', 0.1)
        sp = params.get('spread', 1.0)
        metric = params.get('metric', 'euclidean')
        init = params.get('init', 'spectral')
        param_str = f"umap_n{n}_md{md}_sp{sp}_met{metric}_init{init}"
    else:
        # Generic method parameters
        param_str = method_name.lower() + "_" + "_".join([f"{k}{v}" for k, v in params.items()])
    
    # Create output directory
    param_dir = os.path.join(base_dir, param_str)
    os.makedirs(param_dir, exist_ok=True)
    
    return param_dir


def create_method_factory():
    """
    Create a factory function for dimensional reduction methods
    
    Returns:
    --------
    get_method : function
        Function that returns a dimensional reduction method instance
    """
    def get_method(method_name):
        """
        Get a dimensional reduction method by name
        
        Parameters:
        -----------
        method_name : str
            Name of the method ('tsne', 'umap')
            
        Returns:
        --------
        method : DimensionalReduction
            Instance of a dimensional reduction method
        """
        method_name = method_name.lower()
        
        if method_name == 'tsne':
            from tsne_wrapper import TSNEReduction
            return TSNEReduction()
        elif method_name == 'umap':
            from umap_wrapper import UMAPReduction
            return UMAPReduction()
        else:
            raise ValueError(f"Unknown method: {method_name}. Available methods: 'tsne', 'umap'")
    
    return get_method