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

                            # --- Get the mutation rate and add to expected rates ---
                            # The get_rate method now includes internal debugging
                            rate = self.trinuc_model.get_rate(ref_base, alt_base, context)
                            expected_rates[impact_type][mutation_key] += rate


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

        # Return the dictionaries containing sums of rates and counts per impact type
        # Convert defaultdicts back to regular dicts for safety
        final_expected_rates = {k: dict(v) for k, v in expected_rates.items()}
        final_expected_counts = {k: dict(v) for k, v in expected_counts.items()}

        return {'rates': final_expected_rates, 'counts': final_expected_counts}

