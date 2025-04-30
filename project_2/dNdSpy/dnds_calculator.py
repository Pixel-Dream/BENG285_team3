# dnds_calculator.py - With Debugging
import numpy as np
from typing import Dict, List, Tuple, Optional, Union, Any
from collections import defaultdict
import pandas as pd
import sys
from data_classes import Mutation

# Define how many debug messages to print per category
DEBUG_LIMIT = 5
debug_counters = defaultdict(int)

class DnDsCalculator:
    """
    Calculator for dN/dS ratios (with debugging prints).
    """

    def __init__(self, trinuc_model):
        """Initialize the dN/dS calculator"""
        print("Debug (Calculator.__init__): Initializing DnDsCalculator...")
        if trinuc_model is None:
             raise ValueError("TrinucleotideModel instance is required.")
        self.trinuc_model = trinuc_model

        # Genetic code (standard) - remains the same
        self.genetic_code = {
            'TTT': 'F', 'TTC': 'F', 'TTA': 'L', 'TTG': 'L', 'CTT': 'L',
            'CTC': 'L', 'CTA': 'L', 'CTG': 'L', 'ATT': 'I', 'ATC': 'I',
            'ATA': 'I', 'ATG': 'M', 'GTT': 'V', 'GTC': 'V', 'GTA': 'V',
            'GTG': 'V', 'TCT': 'S', 'TCC': 'S', 'TCA': 'S', 'TCG': 'S',
            'CCT': 'P', 'CCC': 'P', 'CCA': 'P', 'CCG': 'P', 'ACT': 'T',
            'ACC': 'T', 'ACA': 'T', 'ACG': 'T', 'GCT': 'A', 'GCC': 'A',
            'GCA': 'A', 'GCG': 'A', 'TAT': 'Y', 'TAC': 'Y', 'TAA': '*',
            'TAG': '*', 'CAT': 'H', 'CAC': 'H', 'CAA': 'Q', 'CAG': 'Q',
            'AAT': 'N', 'AAC': 'N', 'AAA': 'K', 'AAG': 'K', 'GAT': 'D',
            'GAC': 'D', 'GAA': 'E', 'GAG': 'E', 'TGT': 'C', 'TGC': 'C',
            'TGA': '*', 'TGG': 'W', 'CGT': 'R', 'CGC': 'R', 'CGA': 'R',
            'CGG': 'R', 'AGT': 'S', 'AGC': 'S', 'AGA': 'R', 'AGG': 'R',
            'GGT': 'G', 'GGC': 'G', 'GGA': 'G', 'GGG': 'G'
        }

        # Cache for mutation impact calculations
        self.impact_cache = {}
        print("Debug (Calculator.__init__): DnDsCalculator initialization complete.")
        sys.stdout.flush()


    def _calculate_site_impacts(self, codon: str) -> Dict[Tuple[int, str], str]:
        """
        Calculate the impact of all possible single base substitutions in a codon.
        Returns a dict mapping (pos_in_codon, alt_base) -> impact_type.
        """
        # Use codon directly as cache key for all its site impacts
        if codon in self.impact_cache:
            return self.impact_cache[codon]

        impacts = {}
        ref_aa = self.genetic_code.get(codon, 'X') # 'X' for unknown codons

        if ref_aa == 'X' or len(codon) != 3: # Cannot determine impact if codon is invalid
             self.impact_cache[codon] = {} # Cache empty result
             return {}

        for pos_in_codon in range(3):
            ref_base = codon[pos_in_codon]
            for alt_base in "ACGT":
                if alt_base != ref_base:
                    # Create mutated codon
                    mut_codon_list = list(codon)
                    mut_codon_list[pos_in_codon] = alt_base
                    mut_codon = ''.join(mut_codon_list)

                    # Get resulting amino acid
                    alt_aa = self.genetic_code.get(mut_codon, 'X')

                    # Determine impact
                    impact_type = 'other' # Default
                    if alt_aa == 'X':
                         impact_type = 'other' # Could be invalid codon result
                    elif alt_aa == ref_aa:
                        impact_type = 'synonymous'
                    elif alt_aa == '*' and ref_aa != '*':
                        impact_type = 'nonsense'
                    elif alt_aa != '*' and ref_aa == '*':
                         impact_type = 'stop_loss' # Explicitly label stop loss
                    elif alt_aa != '*' and ref_aa != '*':
                         impact_type = 'missense'
                    # else: impact remains 'other'

                    impacts[(pos_in_codon, alt_base)] = impact_type

        # Cache result for the entire codon
        self.impact_cache[codon] = impacts
        return impacts


    def calculate_expected_ns_ratio(self, sequence: str) -> Dict[str, Dict[str, Dict[str, float]]]:
        """
        Calculate the expected numbers (counts) and summed rates of synonymous and non-synonymous
        mutations based on sequence context and the trinucleotide model.
        Returns: {'rates': {'synonymous': {mut_type: rate_sum}, ...},
                  'counts': {'synonymous': {mut_type: count_sum}, ...}}
        """
        # Reset debug counters for this calculation
        global debug_counters
        # Keep some counters persistent if needed across calls, otherwise reset here
        # debug_counters = defaultdict(int)

        if debug_counters['exp_calc_calls'] < DEBUG_LIMIT:
             print(f"Debug (Calculator.ExpRatio): Calculating expected rates and counts for sequence len={len(sequence)}...")
             debug_counters['exp_calc_calls'] += 1
        sys.stdout.flush()


        # --- Input Validation and Cleaning ---
        if not isinstance(sequence, str) or len(sequence) == 0:
             print("Debug (Calculator.ExpRatio): ERROR - Input sequence is empty or not a string.")
             return {'rates': {'synonymous': {}, 'missense': {}, 'nonsense': {}},
                     'counts': {'synonymous': {}, 'missense': {}, 'nonsense': {}}} # Return empty structure

        # Ensure uppercase and replace non-ACGT with N
        clean_seq = ''.join(n if n in "ACGT" else "N" for n in sequence.upper())
        original_len = len(sequence)
        clean_len = len(clean_seq)
        if clean_len != original_len or clean_seq != sequence.upper():
             if debug_counters['exp_seq_cleaned'] < DEBUG_LIMIT:
                  print(f"  Debug (Calculator.ExpRatio): Sequence cleaned (len {original_len}->{clean_len}).")
                  debug_counters['exp_seq_cleaned'] += 1


        # Ensure sequence length is a multiple of 3 for codon processing
        remainder = clean_len % 3
        if remainder > 0:
            if debug_counters['exp_seq_trimmed'] < DEBUG_LIMIT:
                 print(f"  Debug (Calculator.ExpRatio): Sequence length {clean_len} not multiple of 3. Trimming last {remainder} bases.")
                 debug_counters['exp_seq_trimmed'] += 1
            clean_seq = clean_seq[:-remainder]
            clean_len = len(clean_seq) # Update length

        if clean_len < 3:
            print(f"Debug (Calculator.ExpRatio): Sequence too short (len={clean_len}) after cleaning/trimming. Cannot calculate.")
            return {'rates': {'synonymous': {}, 'missense': {}, 'nonsense': {}},
                     'counts': {'synonymous': {}, 'missense': {}, 'nonsense': {}}} # Return empty structure


        # --- Initialize Expected Rates and Counts ---
        # Store sums of rates per mutation type (e.g., C>T) within each impact class
        expected_rates = {
            'synonymous': defaultdict(float),
            'missense': defaultdict(float),
            'nonsense': defaultdict(float),
            'stop_loss': defaultdict(float), # Track stop loss separately
            'other': defaultdict(float)      # Track other impacts
        }
        # Store counts of potential mutations per mutation type within each impact class
        expected_counts = {
            'synonymous': defaultdict(float), # Use float for consistency, though counts are integers
            'missense': defaultdict(float),
            'nonsense': defaultdict(float),
            'stop_loss': defaultdict(float),
            'other': defaultdict(float)
        }
        codons_processed = 0
        sites_processed = 0
        sites_skipped_n_in_codon = 0
        sites_skipped_n_in_context = 0
        mutations_considered = 0
        mutations_with_zero_rate = 0


        # --- Iterate through Codons ---
        for i in range(0, clean_len - 2, 3): # Iterate by codon start position
            codon = clean_seq[i : i+3]
            codons_processed += 1

            # Skip codons containing 'N'
            if 'N' in codon:
                 sites_skipped_n_in_codon += 3 # Skip all 3 sites in this codon
                 if debug_counters['exp_skip_n_codon'] < DEBUG_LIMIT:
                     print(f"  Debug (Calculator.ExpRatio): Skipping codon '{codon}' at pos {i} due to 'N'.")
                     debug_counters['exp_skip_n_codon'] += 1
                 continue

            # Pre-calculate impacts for all possible SNVs in this codon
            site_impacts = self._calculate_site_impacts(codon) # Returns {(pos, alt): impact}

            # Iterate through each position within the codon
            for pos_in_codon in range(3):
                sites_processed += 1
                ref_base = codon[pos_in_codon]
                abs_pos = i + pos_in_codon # Absolute position in clean_seq

                # Determine trinucleotide context centered on the current position
                # Handle sequence edges carefully
                if abs_pos == 0: # Beginning of sequence
                    context = "N" + clean_seq[abs_pos : abs_pos+2]
                elif abs_pos == clean_len - 1: # End of sequence
                    context = clean_seq[abs_pos-1 : abs_pos+1] + "N"
                else: # Middle of sequence
                    context = clean_seq[abs_pos-1 : abs_pos+2]

                # Skip site if context contains 'N'
                if 'N' in context:
                     sites_skipped_n_in_context += 1
                     if debug_counters['exp_skip_n_context'] < DEBUG_LIMIT:
                         print(f"  Debug (Calculator.ExpRatio): Skipping site at pos {abs_pos} (codon {codon}) due to 'N' in context '{context}'.")
                         debug_counters['exp_skip_n_context'] += 1
                     continue


                # Consider all possible alternative bases
                for alt_base in "ACGT":
                    if alt_base != ref_base:
                        mutations_considered += 1
                        mutation_key = f"{ref_base}>{alt_base}" # e.g., C>T

                        # Get the pre-calculated impact for this specific substitution
                        impact_type = site_impacts.get((pos_in_codon, alt_base))

                        if impact_type: # Should always be found if codon was valid
                            # --- Count the potential mutation ---
                            if impact_type in expected_counts:
                                expected_counts[impact_type][mutation_key] += 1.0
                            else:
                                expected_counts['other'][mutation_key] += 1.0
                                # Log if needed, but impact_type should be known

                            # --- Get the mutation rate and add to expected rates ---
                            # The get_rate method now includes internal debugging
                            rate = self.trinuc_model.get_rate(ref_base, alt_base, context)

                            if rate > 0:
                                # Add the rate to the corresponding impact and mutation type
                                if impact_type in expected_rates:
                                     expected_rates[impact_type][mutation_key] += rate
                                else:
                                     # This case should ideally not happen with current impacts
                                     expected_rates['other'][mutation_key] += rate
                                     if debug_counters['exp_unknown_impact'] < DEBUG_LIMIT:
                                          print(f"  Debug (Calculator.ExpRatio): Unknown impact type '{impact_type}' encountered for {mutation_key} in {codon} ({context}). Added to 'other' rates.")
                                          debug_counters['exp_unknown_impact'] += 1

                            else: # Rate is zero or negative (error?)
                                 mutations_with_zero_rate += 1
                                 # Optional: Debug zero rates (might be very verbose)
                                 # if debug_counters['exp_zero_rate'] < DEBUG_LIMIT:
                                 #      print(f"  Debug (Calculator.ExpRatio): Zero rate returned for {mutation_key} in context {context} (codon {codon}).")
                                 #      debug_counters['exp_zero_rate'] += 1
                        else:
                             # This indicates an issue with _calculate_site_impacts or the codon
                             if debug_counters['exp_impact_not_found'] < DEBUG_LIMIT:
                                  print(f"  Debug (Calculator.ExpRatio): CRITICAL - Impact not found for pos={pos_in_codon}, alt={alt_base} in codon '{codon}'. This shouldn't happen.")
                                  debug_counters['exp_impact_not_found'] += 1


        # --- Summary and Fallback ---
        total_expected_rate_sum = sum(sum(rates.values()) for rates in expected_rates.values())
        total_expected_count_sum = sum(sum(counts.values()) for counts in expected_counts.values())

        if debug_counters['exp_calc_summary'] < DEBUG_LIMIT:
             print(f"Debug (Calculator.ExpRatio): --- Calculation Summary ---")
             print(f"  Codons processed: {codons_processed}")
             print(f"  Sites processed: {sites_processed}")
             print(f"  Sites skipped (N in codon): {sites_skipped_n_in_codon}")
             print(f"  Sites skipped (N in context): {sites_skipped_n_in_context}")
             print(f"  Mutations considered: {mutations_considered}")
             print(f"  Mutations with zero rate: {mutations_with_zero_rate}")
             print(f"  Total expected rate sum: {total_expected_rate_sum:.4f}")
             print(f"  Total expected count sum: {total_expected_count_sum:.1f}")
             print("  Rates Sums per impact:")
             for impact, rates_dict in expected_rates.items():
                  print(f"    - {impact}: {sum(rates_dict.values()):.4f}")
             print("  Counts per impact:")
             for impact, counts_dict in expected_counts.items():
                  print(f"    - {impact}: {sum(counts_dict.values()):.1f}")
             print(f"Debug (Calculator.ExpRatio): --------------------------")
             debug_counters['exp_calc_summary'] += 1
             sys.stdout.flush()


        # Check for the "all zeros" case from the trinucleotide model for rates
        if total_expected_rate_sum == 0 and mutations_considered > 0:
            # This suggests the trinuc model itself is returning zeros, despite valid sites/mutations
            print("Warning (Calculator.ExpRatio): Trinucleotide model returned all zeros. Expected rates sum is 0.")
            print("                             Rates will be zero. Counts are based on potential mutations.")
            # No fallback calculation needed here anymore, as counts are already calculated.
            # The expected_rates dictionary correctly reflects the zero rates.
            sys.stdout.flush()

        # Return the dictionaries containing sums of rates and counts per impact type
        # Convert defaultdicts back to regular dicts for safety
        final_expected_rates = {k: dict(v) for k, v in expected_rates.items()}
        final_expected_counts = {k: dict(v) for k, v in expected_counts.items()}

        return {'rates': final_expected_rates, 'counts': final_expected_counts}


    def calculate_dnds(self, observed_mutations: List[Mutation],
                      expected_output: Dict[str, Dict[str, Dict[str, float]]]) -> Dict[str, float]:
        """
        Calculate dN/dS ratios for a gene based on observed and expected rates/counts.
        Expected output should be the dictionary returned by calculate_expected_ns_ratio.
        (Debugging added)
        """
        if debug_counters['dnds_calc_calls'] < DEBUG_LIMIT:
             print(f"Debug (Calculator.dNdS): Calculating dN/dS for gene...")
             debug_counters['dnds_calc_calls'] += 1
        sys.stdout.flush()

        # --- Extract Expected Rates ---
        # Check if the input dictionary has the expected structure
        if not isinstance(expected_output, dict) or 'rates' not in expected_output:
             print("Error (Calculator.dNdS): Invalid format for expected_output. Expected {'rates': {...}, 'counts': {...}}.")
             # Return NaN for all results or raise an error
             results_nan = {f'dnds_{imp}': np.nan for imp in ['missense', 'nonsense', 'splice_site', 'global']}
             results_nan.update({'n_syn': 0, 'n_mis': 0, 'n_non': 0, 'n_spl': 0,
                                 'exp_syn': np.nan, 'exp_mis': np.nan, 'exp_non': np.nan, 'exp_spl': np.nan,
                                 'exp_cts_syn': np.nan, 'exp_cts_mis': np.nan, 'exp_cts_non': np.nan, 'exp_cts_spl': np.nan}) # Add counts
             return results_nan

        expected_rates = expected_output.get('rates', {}) # Default to empty dict if key missing
        expected_counts_data = expected_output.get('counts', {}) # Get counts too

        # --- Count Observed Mutations by Impact ---
        observed_counts = defaultdict(int)
        for mutation in observed_mutations:
            if mutation.is_synonymous(): observed_counts['synonymous'] += 1
            elif mutation.is_missense(): observed_counts['missense'] += 1
            elif mutation.is_nonsense(): observed_counts['nonsense'] += 1
            elif mutation.is_splice_site(): observed_counts['splice_site'] += 1 # Include splice if relevant
            # Add other categories if needed

        # --- Sum Expected Rates by Impact ---
        # The input 'expected_rates' already contains the summed rates per mutation type
        # We need to sum these rates across mutation types for each impact class
        expected_rate_sums = defaultdict(float)
        if expected_rates: # Check if dict is not empty
             for impact_type, rates_dict in expected_rates.items():
                  expected_rate_sums[impact_type] = sum(rates_dict.values())
        else:
             print("Debug (Calculator.dNdS): WARNING - Expected rates dictionary is empty or missing.")

        # --- Sum Expected Counts by Impact --- (Optional: For reporting/debugging)
        expected_count_sums = defaultdict(float)
        if expected_counts_data:
            for impact_type, counts_dict in expected_counts_data.items():
                expected_count_sums[impact_type] = sum(counts_dict.values())


        print(f"Debug (Calculator.dNdS): --- dN/dS Calculation Inputs ---")
        print(f"  Observed: Syn={observed_counts['synonymous']}, Mis={observed_counts['missense']}, Non={observed_counts['nonsense']}, Spl={observed_counts['splice_site']}")
        print(f"  Expected Rate Sums: Syn={expected_rate_sums['synonymous']:.4f}, Mis={expected_rate_sums['missense']:.4f}, Non={expected_rate_sums['nonsense']:.4f}, Spl={expected_rate_sums['splice_site']:.4f}")
        print(f"  Expected Count Sums: Syn={expected_count_sums['synonymous']:.1f}, Mis={expected_count_sums['missense']:.1f}, Non={expected_count_sums['nonsense']:.1f}, Spl={expected_count_sums['splice_site']:.1f}") # Print counts


        # --- Calculate dN/dS Ratios ---
        # dN/dS fundamentally relies on RATES (observed count / expected rate)
        dnds_results = {}
        n_syn_obs = observed_counts['synonymous']
        e_syn_exp = expected_rate_sums['synonymous'] # Use EXPECTED RATE sum

        # Check if synonymous baseline can be calculated
        if n_syn_obs > 0 and e_syn_exp > 0:
            syn_rate = n_syn_obs / e_syn_exp # Observed synonymous rate
            print(f"Debug (Calculator.dNdS): Synonymous Rate (Obs/ExpRate) = {syn_rate:.4f}")

            # Calculate dN/dS for different non-synonymous types
            for impact_type in ['missense', 'nonsense', 'splice_site']: # Add others if needed
                 n_obs = observed_counts[impact_type]
                 e_exp = expected_rate_sums[impact_type] # Use EXPECTED RATE sum

                 if e_exp > 0:
                      non_syn_rate = n_obs / e_exp
                      dnds_results[f'dnds_{impact_type}'] = non_syn_rate / syn_rate
                      print(f"  Debug (Calculator.dNdS): {impact_type.capitalize()} Rate={non_syn_rate:.4f}, dN/dS={dnds_results[f'dnds_{impact_type}']:.4f}")
                 else:
                      dnds_results[f'dnds_{impact_type}'] = np.nan # Assign NaN if expected rate is zero
                      print(f"  Debug (Calculator.dNdS): {impact_type.capitalize()} Rate=N/A (ExpRate=0), dN/dS=NaN")


            # Calculate Global dN/dS (all non-synonymous combined)
            n_nonsyn_obs = observed_counts['missense'] + observed_counts['nonsense'] + observed_counts['splice_site']
            e_nonsyn_exp = expected_rate_sums['missense'] + expected_rate_sums['nonsense'] + expected_rate_sums['splice_site'] # Use EXPECTED RATE sums

            if e_nonsyn_exp > 0:
                 global_nonsyn_rate = n_nonsyn_obs / e_nonsyn_exp
                 dnds_results['dnds_global'] = global_nonsyn_rate / syn_rate
                 print(f"  Debug (Calculator.dNdS): Global NonSyn Rate={global_nonsyn_rate:.4f}, dN/dS={dnds_results['dnds_global']:.4f}")
            else:
                 dnds_results['dnds_global'] = np.nan
                 print(f"  Debug (Calculator.dNdS): Global NonSyn Rate=N/A (ExpRate=0), dN/dS=NaN")

        else:
            # Cannot calculate any dN/dS if baseline is zero
            print("Debug (Calculator.dNdS): WARNING - Cannot calculate dN/dS (ObsSyn=0 or ExpRateSyn=0). Setting all to NaN.")
            for impact_type in ['missense', 'nonsense', 'splice_site', 'global']:
                 dnds_results[f'dnds_{impact_type}'] = np.nan


        # Add raw counts, expected rate sums, AND expected count sums to the results
        dnds_results.update({
            'n_syn': observed_counts['synonymous'],
            'n_mis': observed_counts['missense'],
            'n_non': observed_counts['nonsense'],
            'n_spl': observed_counts['splice_site'],
            'exp_rate_syn': expected_rate_sums['synonymous'],
            'exp_rate_mis': expected_rate_sums['missense'],
            'exp_rate_non': expected_rate_sums['nonsense'],
            'exp_rate_spl': expected_rate_sums['splice_site'],
            'exp_cts_syn': expected_count_sums['synonymous'], # Add expected counts
            'exp_cts_mis': expected_count_sums['missense'],
            'exp_cts_non': expected_count_sums['nonsense'],
            'exp_cts_spl': expected_count_sums['splice_site']
        })

        print("Debug (Calculator.dNdS): Finished dN/dS calculation for gene.")
        sys.stdout.flush()
        return dnds_results # Return the dictionary with dN/dS values and counts

