import numpy as np
import pandas as pd
import time
from umap.umap_ import UMAP
from dimension_reduction import DimensionalReduction

class UMAPReduction(DimensionalReduction):
    """UMAP implementation of dimensional reduction."""
    
    def fit_transform(self, data, n_neighbors=15, min_dist=0.1, n_components=2,
                     metric='euclidean', spread=1.0, n_epochs=None, 
                     learning_rate=1.0, init='spectral', random_state=42, 
                     low_memory=False, verbose=False, **kwargs):
        """
        Run UMAP on preprocessed data
        
        Parameters:
        -----------
        data : DataFrame
            Preprocessed expression data with samples as rows
        n_neighbors : int
            Number of neighbors to consider (5-50)
        min_dist : float
            Minimum distance between points in embedding (0.0-1.0)
        n_components : int
            Number of dimensions for embedding
        metric : str
            Distance metric ('euclidean', 'correlation', 'cosine', etc.)
        spread : float
            Spread of the embedding
        n_epochs : int or None
            Number of epochs for embedding
        learning_rate : float
            Initial learning rate
        init : str
            Initialization method ('spectral', 'random')
        random_state : int
            Random seed for reproducibility
        low_memory : bool
            Use memory-efficient algorithm variant
        verbose : bool
            Verbosity level
            
        Returns:
        --------
        umap_result : DataFrame
            DataFrame with UMAP coordinates
        """
        print(f"Running UMAP with parameters: n_neighbors={n_neighbors}, "
              f"min_dist={min_dist}, spread={spread}, metric={metric}, "
              f"random_state={random_state}")
        
        start_time = time.time()
        
        reducer = UMAP(
            n_neighbors=n_neighbors,
            min_dist=min_dist,
            n_components=n_components,
            metric=metric,
            spread=spread,
            n_epochs=n_epochs,
            learning_rate=learning_rate,
            init=init,
            random_state=random_state,
            low_memory=low_memory,
            verbose=verbose
        )
        
        umap_embedding = reducer.fit_transform(data)
        
        elapsed_time = time.time() - start_time
        print(f"UMAP completed in {elapsed_time:.2f} seconds")
        
        # Create a DataFrame with UMAP results
        umap_result = pd.DataFrame(
            umap_embedding,
            index=data.index,
            columns=['UMAP1', 'UMAP2'] if n_components == 2 else [f'UMAP{i+1}' for i in range(n_components)]
        )
        
        return umap_result
    
    def get_method_name(self):
        """Return the name of the method."""
        return "UMAP"
    
    def get_default_params(self):
        """Return default parameters for the method."""
        return {
            'n_neighbors': 15,
            'min_dist': 0.1,
            'n_components': 2,
            'metric': 'euclidean',
            'spread': 1.0,
            'n_epochs': None,
            'learning_rate': 1.0,
            'init': 'spectral',
            'verbose': False
        }
    
    def get_param_grid(self):
        """Return a reasonable parameter grid for exploration."""
        return {
            'n_neighbors': [5, 15, 30],  # Added higher value
            'min_dist': [0.0, 0.1, 0.5],  # Added higher value
            'spread': [1.0, 1.5, 2.0],  # Added options
            'metric': ['euclidean', 'cosine'],  # Added alternatives
            'init': ['spectral', 'random', 'pca'],  # Added random initialization
            'n_epochs': [200, 500],  # Added epoch options
            'learning_rate': [0.5, 1.0],  # Added learning rate options
            'random_state': [123]  # Added for reproducibility
        }