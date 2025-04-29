# negative_binomial_model.py
import numpy as np
import pandas as pd
import statsmodels.api as sm
from typing import Dict, List, Optional, Union, Any, Tuple
from scipy.stats import nbinom

class NegativeBinomialModel:
    """
    Negative Binomial regression model for background mutation rates.
    
    This model accounts for the overdispersion of mutation rates across genes,
    which is a key feature of cancer genomics data. It incorporates covariates
    like gene expression, replication timing, and chromatin state.
    """
    
    def __init__(self):
        self.theta = None  # Dispersion parameter
        self.coefficients = None  # Regression coefficients
        self.model = None  # Fitted statsmodels model
    
    def fit(self, genes: Dict[str, Any], covariates_df: pd.DataFrame = None) -> 'NegativeBinomialModel':
        """
        Fit a negative binomial regression model
        
        Parameters:
        -----------
        genes: Dict[str, Any]
            Dictionary of Gene objects
        covariates_df: pd.DataFrame
            DataFrame of covariates with gene names as index
            
        Returns:
        --------
        self: NegativeBinomialModel
            Returns self for method chaining
        """
        # Prepare data
        gene_data = []
        for gene_name, gene in genes.items():
            if gene_name in covariates_df.index and gene.coding_length > 0:
                # Extract gene data
                row = {
                    'gene_name': gene_name,
                    'n_syn': gene.count_synonymous(),
                    'expected_syn': gene.coding_length  # Simplified - will be refined
                }
                
                # Add covariates
                for col in covariates_df.columns:
                    row[col] = covariates_df.loc[gene_name, col]
                
                gene_data.append(row)
        
        if not gene_data:
            raise ValueError("No valid genes with covariates for model fitting")
        
        # Convert to DataFrame
        df = pd.DataFrame(gene_data)
        
        # Use the log of expected synonymous mutations as offset
        df['log_exp_syn'] = np.log(df['expected_syn'])
        
        # Exclude genes with zero expected synonymous mutations
        df = df[df['expected_syn'] > 0]
        
        # Prepare formula for statsmodels
        if covariates_df is not None and covariates_df.shape[1] > 0:
            formula = 'n_syn ~ ' + ' + '.join(covariates_df.columns)
        else:
            formula = 'n_syn ~ 1'  # Intercept only
            
        # Add offset term
        formula += ' + offset(log_exp_syn)'
        
        # Fit negative binomial model
        try:
            model = sm.formula.negativebinomial.NegativeBinomial(formula, df).fit(disp=0)
            self.model = model
            self.coefficients = model.params
            self.theta = model.alpha  # Dispersion parameter
        except Exception as e:
            # Fallback to simpler model if fitting fails
            print(f"Negative binomial fitting failed: {e}")
            print("Falling back to simpler model")
            formula = 'n_syn ~ 1 + offset(log_exp_syn)'
            model = sm.formula.negativebinomial.NegativeBinomial(formula, df).fit(disp=0)
            self.model = model
            self.coefficients = model.params
            self.theta = model.alpha
        
        return self
    
    def predict_mutation_rate(self, gene_name: str, covariates: Dict[str, float], 
                             coding_length: int) -> float:
        """
        Predict the mutation rate for a gene
        
        Parameters:
        -----------
        gene_name: str
            Name of the gene
        covariates: Dict[str, float]
            Dictionary of covariate values
        coding_length: int
            Length of coding sequence
            
        Returns:
        --------
        rate: float
            Predicted mutation rate
        """
        if self.model is None:
            raise ValueError("Model not fitted yet")
        
        # Create a DataFrame with covariates for prediction
        row = pd.DataFrame([covariates], index=[gene_name])
        row['log_exp_syn'] = np.log(coding_length)
        
        # Predict using the model
        predicted_rate = self.model.predict(row).iloc[0]
        return predicted_rate
    
    def get_gamma_params(self, gene_name: str, covariates: Dict[str, float], 
                        coding_length: int) -> Tuple[float, float]:
        """
        Get the parameters of the Gamma distribution for mutation rate
        
        Parameters:
        -----------
        gene_name: str
            Name of the gene
        covariates: Dict[str, float]
            Dictionary of covariate values
        coding_length: int
            Length of coding sequence
            
        Returns:
        --------
        shape: float
            Shape parameter of the Gamma distribution
        scale: float
            Scale parameter of the Gamma distribution
        """
        if self.theta is None:
            raise ValueError("Model not fitted yet")
        
        # Predict mean
        mu = self.predict_mutation_rate(gene_name, covariates, coding_length)
        
        # Calculate Gamma parameters
        shape = self.theta
        scale = mu / self.theta
        
        return shape, scale
    
    def save(self, filename: str) -> None:
        """Save the model to a file"""
        import pickle
        # Don't save the statsmodels object directly
        to_save = {
            'theta': self.theta,
            'coefficients': self.coefficients
        }
        with open(filename, 'wb') as f:
            pickle.dump(to_save, f)
    
    @classmethod
    def load(cls, filename: str) -> 'NegativeBinomialModel':
        """Load a model from a file"""
        import pickle
        model = cls()
        with open(filename, 'rb') as f:
            saved = pickle.load(f)
            model.theta = saved['theta']
            model.coefficients = saved['coefficients']
        return model
    