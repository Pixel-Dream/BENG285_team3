# selection_tester.py
import numpy as np
import pandas as pd
from scipy import stats
from typing import Dict, List, Tuple, Optional, Union, Any
from statsmodels.stats.multitest import multipletests

class SelectionTester:
    """
    Statistical testing module for detecting selection in cancer genomes.
    
    This class implements likelihood ratio tests for detecting genes under
    positive or negative selection, with appropriate multiple testing correction.
    """
    
    def __init__(self):
        pass
    
    def likelihood_ratio_test(self, observed: Dict[str, int], expected: Dict[str, float], 
                            null_model: str = 'neutral') -> Tuple[float, float]:
        """
        Perform likelihood ratio test for selection
        
        Parameters:
        -----------
        observed: Dict[str, int]
            Observed mutation counts by type
        expected: Dict[str, float]
            Expected mutation counts by type
        null_model: str
            Null model to test against ('neutral', 'positive', 'negative')
            
        Returns:
        --------
        likelihood_ratio: float
            Likelihood ratio statistic
        p_value: float
            P-value from chi-squared distribution
        """
        # Extract counts
        n_syn = observed.get('synonymous', 0)
        n_mis = observed.get('missense', 0)
        n_non = observed.get('nonsense', 0)
        n_spl = observed.get('splice_site', 0)
        
        e_syn = expected.get('synonymous', 0.0)
        e_mis = expected.get('missense', 0.0)
        e_non = expected.get('nonsense', 0.0)
        e_spl = expected.get('splice_site', 0.0)
        
        # Skip genes with insufficient mutations
        if n_syn + n_mis + n_non + n_spl == 0 or e_syn + e_mis + e_non + e_spl == 0:
            return 0.0, 1.0
        
        # Calculate total non-synonymous
        n_nonsyn = n_mis + n_non + n_spl
        e_nonsyn = e_mis + e_non + e_spl
        
        # Calculate null and alternative log-likelihoods
        if null_model == 'neutral':
            # Null: dN/dS = 1 (neutral evolution)
            # Alternative: dN/dS != 1 (selection)
            
            # Null model (dN/dS = 1)
            mu_null = (n_syn + n_nonsyn) / (e_syn + e_nonsyn)
            ll_null = stats.poisson.logpmf(n_syn, mu_null * e_syn) + stats.poisson.logpmf(n_nonsyn, mu_null * e_nonsyn)
            
            # Alternative model (separate rates)
            mu_syn = n_syn / e_syn if e_syn > 0 else 0
            mu_nonsyn = n_nonsyn / e_nonsyn if e_nonsyn > 0 else 0
            ll_alt = stats.poisson.logpmf(n_syn, mu_syn * e_syn) + stats.poisson.logpmf(n_nonsyn, mu_nonsyn * e_nonsyn)
            
            # Likelihood ratio test
            lr = -2 * (ll_null - ll_alt)
            p_value = stats.chi2.sf(lr, df=1)  # 1 degree of freedom
            
        elif null_model == 'negative':
            # Null: dN/dS <= 1 (negative or neutral)
            # Alternative: dN/dS > 1 (positive selection)
            
            if n_nonsyn / e_nonsyn <= n_syn / e_syn:
                # If observed dN/dS <= 1, null model is best fit
                return 0.0, 1.0
            
            # Null model (dN/dS = 1)
            mu_null = (n_syn + n_nonsyn) / (e_syn + e_nonsyn)
            ll_null = stats.poisson.logpmf(n_syn, mu_null * e_syn) + stats.poisson.logpmf(n_nonsyn, mu_null * e_nonsyn)
            
            # Alternative model (dN/dS > 1)
            mu_syn = n_syn / e_syn if e_syn > 0 else 0
            mu_nonsyn = n_nonsyn / e_nonsyn if e_nonsyn > 0 else 0
            ll_alt = stats.poisson.logpmf(n_syn, mu_syn * e_syn) + stats.poisson.logpmf(n_nonsyn, mu_nonsyn * e_nonsyn)
            
            # One-sided test
            lr = -2 * (ll_null - ll_alt)
            p_value = 0.5 * stats.chi2.sf(lr, df=1)  # One-sided test
            
        elif null_model == 'positive':
            # Null: dN/dS >= 1 (positive or neutral)
            # Alternative: dN/dS < 1 (negative selection)
            
            if n_nonsyn / e_nonsyn >= n_syn / e_syn:
                # If observed dN/dS >= 1, null model is best fit
                return 0.0, 1.0
            
            # Null model (dN/dS = 1)
            mu_null = (n_syn + n_nonsyn) / (e_syn + e_nonsyn)
            ll_null = stats.poisson.logpmf(n_syn, mu_null * e_syn) + stats.poisson.logpmf(n_nonsyn, mu_null * e_nonsyn)
            
            # Alternative model (dN/dS < 1)
            mu_syn = n_syn / e_syn if e_syn > 0 else 0
            mu_nonsyn = n_nonsyn / e_nonsyn if e_nonsyn > 0 else 0
            ll_alt = stats.poisson.logpmf(n_syn, mu_syn * e_syn) + stats.poisson.logpmf(n_nonsyn, mu_nonsyn * e_nonsyn)
            
            # One-sided test
            lr = -2 * (ll_null - ll_alt)
            p_value = 0.5 * stats.chi2.sf(lr, df=1)  # One-sided test
        
        else:
            raise ValueError(f"Unknown null model: {null_model}")
        
        return lr, p_value
    
    def test_selection(self, genes: Dict[str, Any], observed_mutations: Dict[str, List], 
                     expected_counts_rates: Dict[str, Dict[str, Dict[str, float]]], 
                     expected_counts_counts: Dict[str, Dict[str, Dict[str, float]]], 
                     null_model: str = 'neutral', fdr_threshold: float = 0.1) -> pd.DataFrame:
        """
        Test for selection across genes
        
        Parameters:
        -----------
        genes: Dict[str, Any]
            Dictionary of Gene objects
        observed_mutations: Dict[str, List]
            Dictionary mapping gene names to lists of observed mutations
        expected_counts_rates: Dict[str, Dict[str, Dict[str, float]]]
            Dictionary mapping gene names to expected mutation rates by type
        expected_counts_counts: Dict[str, Dict[str, Dict[str, float]]]
            Dictionary mapping gene names to expected mutation counts by type
        null_model: str
            Null model to test against ('neutral', 'positive', 'negative')
        fdr_threshold: float
            False discovery rate threshold for significance
            
        Returns:
        --------
        results: pd.DataFrame
            DataFrame with test results for each gene
        """
        results = []
        
        for gene_name, gene in genes.items():
            if gene_name not in expected_counts_rates:
                continue
                
            # Count observed mutations by type
            observed = {
                'synonymous': 0,
                'missense': 0,
                'nonsense': 0,
                'splice_site': 0
            }
            
            if gene_name in observed_mutations:
                for mutation in observed_mutations[gene_name]:
                    if mutation.is_synonymous():
                        observed['synonymous'] += 1
                    elif mutation.is_missense():
                        observed['missense'] += 1
                    elif mutation.is_nonsense():
                        observed['nonsense'] += 1
                    elif mutation.is_splice_site():
                        observed['splice_site'] += 1
            
            # Calculate total expected counts
            expected = {
                'synonymous': sum(expected_counts_rates[gene_name]['synonymous'].values()),
                'missense': sum(expected_counts_rates[gene_name]['missense'].values()),
                'nonsense': sum(expected_counts_rates[gene_name]['nonsense'].values()),
                'splice_site': sum(expected_counts_rates[gene_name].get('splice_site', {}).values())
            }

            expected_counts = {
                'synonymous': sum(expected_counts_counts[gene_name]['synonymous'].values()),
                'missense': sum(expected_counts_counts[gene_name]['missense'].values()),
                'nonsense': sum(expected_counts_counts[gene_name]['nonsense'].values()),
                'splice_site': sum(expected_counts_counts[gene_name].get('splice_site', {}).values())
            }
            
            # Perform likelihood ratio test
            lr, p_value = self.likelihood_ratio_test(observed, expected, null_model)
            
            # Calculate dN/dS ratios
            total_nonsyn = observed['missense'] + observed['nonsense'] + observed['splice_site']
            total_exp_nonsyn = expected['missense'] + expected['nonsense'] + expected['splice_site']
            
            if observed['synonymous'] > 0 and expected['synonymous'] > 0:
                dnds_global = (total_nonsyn / total_exp_nonsyn) / (observed['synonymous'] / expected['synonymous'])
                dnds_mis = (observed['missense'] / expected['missense']) / (observed['synonymous'] / expected['synonymous']) if expected['missense'] > 0 else float('nan')
                dnds_non = (observed['nonsense'] / expected['nonsense']) / (observed['synonymous'] / expected['synonymous']) if expected['nonsense'] > 0 else float('nan')
            else:
                dnds_global = dnds_mis = dnds_non = float('nan')
            
            # Store results
            results.append({
                'gene_name': gene_name,
                'n_syn': observed['synonymous'],
                'n_mis': observed['missense'],
                'n_non': observed['nonsense'],
                'n_spl': observed['splice_site'],
                'exp_syn': expected['synonymous'],
                'exp_mis': expected['missense'],
                'exp_non': expected['nonsense'],
                'exp_spl': expected['splice_site'],
                'exp_syn_counts': expected_counts['synonymous'],
                'exp_mis_counts': expected_counts['missense'],
                'exp_non_counts': expected_counts['nonsense'],
                'exp_spl_counts': expected_counts['splice_site'],
                'dnds_global': dnds_global,
                'dnds_mis': dnds_mis,
                'dnds_non': dnds_non,
                'lr': lr,
                'p_value': p_value
            })
        
        # Convert to DataFrame
        results_df = pd.DataFrame(results)
        
        # Apply multiple testing correction
        if not results_df.empty and 'p_value' in results_df.columns:
            _, q_values, _, _ = multipletests(results_df['p_value'], method='fdr_bh')
            results_df['q_value'] = q_values
            
            # Flag significant genes
            results_df['significant'] = results_df['q_value'] < fdr_threshold
            
            # Sort by p-value
            results_df = results_df.sort_values('p_value')
            
        return results_df
    