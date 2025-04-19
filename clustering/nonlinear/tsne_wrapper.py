import numpy as np
import pandas as pd
from sklearn.manifold import TSNE
import time
from dimension_reduction import DimensionalReduction

class TSNEReduction(DimensionalReduction):
    """t-SNE implementation of dimensional reduction."""
    
    def fit_transform(self, data, perplexity=30, learning_rate=200, n_iter=1000, 
                     early_exaggeration=12, init='pca', metric='euclidean', 
                     random_state=42, verbose=1, **kwargs):
        """
        Run t-SNE on preprocessed data
        
        Parameters:
        -----------
        data : DataFrame
            Preprocessed expression data with samples as rows
        perplexity : float
            Perplexity parameter for t-SNE (5-50)
        learning_rate : float
            Learning rate for t-SNE (100-1000)
        n_iter : int
            Number of iterations (250-1000+)
        early_exaggeration : float
            Early exaggeration factor (4-12)
        init : str
            Initialization method ('random', 'pca', or ndarray)
        metric : str
            Distance metric ('euclidean', 'manhattan', 'cosine', etc.)
        random_state : int
            Random seed for reproducibility
        verbose : int
            Verbosity level
            
        Returns:
        --------
        tsne_result : DataFrame
            DataFrame with t-SNE coordinates
        """
        print(f"Running t-SNE with parameters: perplexity={perplexity}, "
              f"learning_rate={learning_rate}, n_iter={n_iter}, "
              f"early_exaggeration={early_exaggeration}, init={init}, metric={metric}, "
              f"random_state={random_state}")
        
        start_time = time.time()
        
        tsne = TSNE(
            n_components=2,
            perplexity=perplexity,
            learning_rate=learning_rate,
            n_iter=n_iter,
            early_exaggeration=early_exaggeration,
            init=init,
            metric=metric,
            random_state=random_state,
            verbose=verbose
        )
        
        tsne_embedding = tsne.fit_transform(data)
        
        elapsed_time = time.time() - start_time
        print(f"t-SNE completed in {elapsed_time:.2f} seconds")
        
        # Create a DataFrame with t-SNE results
        tsne_result = pd.DataFrame(
            tsne_embedding,
            index=data.index,
            columns=['tSNE1', 'tSNE2']
        )
        
        return tsne_result
    
    def get_method_name(self):
        """Return the name of the method."""
        return "t-SNE"
    
    def get_default_params(self):
        """Return default parameters for the method."""
        return {
            'perplexity': 30,
            'learning_rate': 200,
            'n_iter': 1000,
            'early_exaggeration': 12,
            'init': 'pca',
            'metric': 'euclidean',
            'verbose': 1
        }
    
    def get_param_grid(self):
        """Return a comprehensive parameter grid for t-SNE exploration."""
        return {
            # Your existing parameters (all good choices)
            'perplexity': [10, 30, 50],
            'learning_rate': [100, 200, 500],
            'n_iter': [500, 1000, 2000],
            'early_exaggeration': [8, 12, 16],
            'init': ['pca', 'random'],
            'metric': ['euclidean', 'cosine', 'correlation']
        }