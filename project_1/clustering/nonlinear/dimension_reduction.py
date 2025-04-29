from abc import ABC, abstractmethod
import pandas as pd

class DimensionalReduction(ABC):
    """Base class for dimensional reduction techniques."""
    
    @abstractmethod
    def fit_transform(self, data, **kwargs):
        """
        Run dimensionality reduction and return embedding.
        
        Parameters:
        -----------
        data : DataFrame
            Preprocessed data with samples as rows
        **kwargs : dict
            Method-specific parameters
        
        Returns:
        --------
        embedding : DataFrame
            DataFrame with embedding coordinates
        """
        pass
    
    @abstractmethod
    def get_method_name(self):
        """
        Return the name of the method.
        
        Returns:
        --------
        method_name : str
            Name of the dimensionality reduction method
        """
        pass
    
    @abstractmethod
    def get_default_params(self):
        """
        Return default parameters for the method.
        
        Returns:
        --------
        params : dict
            Dictionary of default parameters
        """
        pass
    
    @abstractmethod
    def get_param_grid(self):
        """
        Return a reasonable parameter grid for exploration.
        
        Returns:
        --------
        param_grid : dict
            Dictionary mapping parameter names to lists of values
        """
        pass