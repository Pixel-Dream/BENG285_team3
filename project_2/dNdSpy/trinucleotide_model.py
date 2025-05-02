# trinucleotide_model.py - Fixed Context Key Logic

import numpy as np
from typing import Dict, List, Tuple, Optional, Union
from collections import defaultdict
import sys # For flushing output
from statsmodels.genmod.generalized_linear_model import GLMResults


# Define how many debug messages to print per category
DEBUG_LIMIT = 10 # Keep debug prints for now
debug_counters = defaultdict(int)

class TrinucleotideModel:
    """
    Model for trinucleotide context-dependent mutation rates.
    (Fixed context key logic + debugging v2)
    """

    def __init__(self):
        # Define substitution types and contexts
        self.sub_types = ['C>A', 'C>G', 'C>T', 'T>A', 'T>C', 'T>G']
        # The context keys EXPECT 'N' in the middle
        self.contexts = [f"{five}N{three}" for five in "ACGT" for three in "ACGT"]

        # Initialize rate matrix (6 substitution types × 16 contexts × 2 strands)
        self.rates = np.ones((6, 16, 2))  # Default to uniform rates

        # Mapping from nucleotide to complementary nucleotide
        self.complement = {
            'A': 'T', 'C': 'G', 'G': 'C', 'T': 'A', 'N': 'N'
        }

        # Create conversion dictionaries
        self._init_conversion_maps()

        # Cache to store results
        self.cache = {}
        print("Debug (TrinucModel): Initialized TrinucleotideModel.")

    def _init_conversion_maps(self):
        """Initialize maps for converting between mutation types and indices"""
        self.sub_to_idx = {sub: i for i, sub in enumerate(self.sub_types)}
        # ctx_to_idx keys are like 'ANA', 'CNC', etc.
        self.ctx_to_idx = {ctx: i for i, ctx in enumerate(self.contexts)}
        print(f"Debug (TrinucModel): ctx_to_idx map initialized with {len(self.ctx_to_idx)} keys (e.g., {list(self.ctx_to_idx.keys())[:5]}...).")


    def _get_complementary_mutation(self, ref: str, alt: str, context: str) -> Tuple[str, str, str]:
        """Get the complementary mutation for a given mutation"""
        comp_ref = self.complement.get(ref, 'N')
        comp_alt = self.complement.get(alt, 'N')
        # Ensure context is treated as string
        comp_context = ''.join([self.complement.get(n, 'N') for n in reversed(str(context))])
        return comp_ref, comp_alt, comp_context

    def _normalize_mutation(self, ref: str, alt: str, context: str) -> Tuple[str, str, str, int]:
        """Normalize a mutation to C>* or T>*"""
        # Ensure context is a string and has length 3 before checking bases
        if not isinstance(context, str) or len(context) != 3:
             return ref, alt, context, -1 # Invalid context format

        # Skip non-standard bases immediately
        if ref not in "ACGT" or alt not in "ACGT" or \
           context[0] not in "ACGTN" or context[2] not in "ACGTN":
             # Strand -1 indicates skip
             return ref, alt, context, -1

        # Ensure we have a standard substitution (C>* or T>*)
        if ref in ['C', 'T']:
            return ref, alt, context, 0 # Strand 0
        else:
            # Convert to complementary mutation
            comp_ref, comp_alt, comp_context = self._get_complementary_mutation(ref, alt, context)
            # Check if complementary is valid before returning
            if comp_ref not in "ACGT" or comp_alt not in "ACGT":
                 return comp_ref, comp_alt, comp_context, -1 # Skip if complement is invalid
            # Also check complementary context format
            if not isinstance(comp_context, str) or len(comp_context) != 3:
                 return comp_ref, comp_alt, comp_context, -1 # Skip if complement context is invalid
            return comp_ref, comp_alt, comp_context, 1 # Strand 1

    def fit(self, mutations: List, ref_genome) -> 'TrinucleotideModel':
        """
        Fit the trinucleotide model to observed mutations (fixed context key logic)
        """
        print("Debug (TrinucModel.fit): Starting model fitting...")
        counts = np.zeros((6, 16, 2))
        processed_mutations = 0
        skipped_context = 0
        skipped_non_snv_acgt = 0
        skipped_normalize = 0
        skipped_sub_key = 0
        skipped_ctx_key = 0
        total_to_process = len(mutations)

        # Reset debug counters for this fit run
        global debug_counters
        debug_counters = defaultdict(int) # Reset counters

        print(f"Debug (TrinucModel.fit): Processing {total_to_process} mutations...")
        sys.stdout.flush() # Ensure prints appear immediately

        # Count observed mutations
        for idx, mutation in enumerate(mutations):
            context = mutation.get_trinucleotide_context(ref_genome)

            if not context or not isinstance(context, str) or len(context) != 3:
                skipped_context += 1
                if debug_counters['fit_skip_context'] < DEBUG_LIMIT:
                    print(f"  >>> Skipping mutation {idx}/{total_to_process} ({mutation}): Invalid context fetched or returned: '{context}' (Type: {type(context)})")
                    debug_counters['fit_skip_context'] += 1
                    sys.stdout.flush()
                continue # Skip this mutation

            # --- Check: Ref/Alt Alleles (SNV, ACGT) ---
            ref = mutation.ref_allele
            alt = mutation.alt_allele
            if not isinstance(ref, str) or not isinstance(alt, str) or \
               len(ref) != 1 or len(alt) != 1 or ref not in "ACGT" or alt not in "ACGT":
                skipped_non_snv_acgt += 1
                if debug_counters['fit_skip_non_snv_acgt'] < DEBUG_LIMIT:
                    print(f"  >>> Skipping mutation {idx}/{total_to_process} ({mutation}): Not SNV or non-ACGT ref/alt. Ref='{ref}', Alt='{alt}'")
                    debug_counters['fit_skip_non_snv_acgt'] += 1
                    sys.stdout.flush()
                continue # Skip this mutation

            # --- Check: Normalization (Handles ambiguous context bases too) ---
            norm_ref, norm_alt, norm_context, strand = self._normalize_mutation(ref, alt, context)
            if strand == -1:
                skipped_normalize += 1
                if debug_counters['fit_skip_normalize'] < DEBUG_LIMIT:
                    print(f"  >>> Skipping mutation {idx}/{total_to_process} ({mutation}): Failed normalization (non-ACGT base/context or invalid complement). Ref='{ref}', Alt='{alt}', Context='{context}' -> NormRef='{norm_ref}', NormAlt='{norm_alt}', NormCtx='{norm_context}'")
                    debug_counters['fit_skip_normalize'] += 1
                    sys.stdout.flush()
                continue # Skip this mutation

            # --- Check: Substitution Key ---
            sub_key = f"{norm_ref}>{norm_alt}"
            if sub_key not in self.sub_to_idx:
                skipped_sub_key += 1
                if debug_counters['fit_skip_sub_key'] < DEBUG_LIMIT:
                    print(f"  >>> Skipping mutation {idx}/{total_to_process} ({mutation}): Invalid substitution key after normalization: '{sub_key}' (NormRef='{norm_ref}', NormAlt='{norm_alt}')")
                    debug_counters['fit_skip_sub_key'] += 1
                    sys.stdout.flush()
                continue # Skip this mutation

            # --- Check: Context Key ---
            # Context key uses flanking bases from normalized context and 'N' in middle
            if not isinstance(norm_context, str) or len(norm_context) != 3:
                # This case might be caught by normalization check, but double check
                skipped_ctx_key += 1
                if debug_counters['fit_skip_ctx_key_len'] < DEBUG_LIMIT:
                    print(f"  >>> Skipping mutation {idx}/{total_to_process} ({mutation}): Invalid normalized context format/length: '{norm_context}'")
                    debug_counters['fit_skip_ctx_key_len'] += 1
                    sys.stdout.flush()
                continue

            ctx_key = f"{norm_context[0]}N{norm_context[2]}"
            ctx_idx = self.ctx_to_idx.get(ctx_key)

            # --- If all checks passed ---
            sub_idx = self.sub_to_idx[sub_key]
            # Ensure indices are valid before incrementing (should be if keys were found)
            if 0 <= sub_idx < counts.shape[0] and 0 <= ctx_idx < counts.shape[1] and 0 <= strand < counts.shape[2]:
                # MAIN LOGIC
                counts[sub_idx, ctx_idx, strand] += 1
                processed_mutations += 1
                # Print first few successful ones
                if debug_counters['fit_success'] < DEBUG_LIMIT:
                    print(f"  --- Processed mutation {idx}/{total_to_process} ({mutation}): Context='{context}', NormRef='{norm_ref}', NormAlt='{norm_alt}', NormCtx='{norm_context}', Strand={strand}, SubKey='{sub_key}'(idx={sub_idx}), CtxKey='{ctx_key}'(idx={ctx_idx})")
                    debug_counters['fit_success'] += 1
                    sys.stdout.flush()

            else:
                 # This would indicate a logic error earlier
                 print(f"  >>> CRITICAL ERROR mutation {idx}/{total_to_process} ({mutation}): Calculated indices out of bounds! sub={sub_idx}, ctx={ctx_idx}, strand={strand}")
                 sys.stdout.flush()

        print(f"\nDebug (TrinucModel.fit): --- Fitting Summary ---")
        print(f"Debug (TrinucModel.fit): Total mutations input: {total_to_process}")
        print(f"Debug (TrinucModel.fit): Successfully processed: {processed_mutations}")
        print(f"Debug (TrinucModel.fit): Skipped - Invalid context fetch/format: {skipped_context}")
        print(f"Debug (TrinucModel.fit): Skipped - Not SNV / Non-ACGT ref/alt: {skipped_non_snv_acgt}")
        print(f"Debug (TrinucModel.fit): Skipped - Failed normalization/Invalid context: {skipped_normalize}")
        print(f"Debug (TrinucModel.fit): Skipped - Invalid substitution key: {skipped_sub_key}")
        print(f"Debug (TrinucModel.fit): Skipped - Invalid context key (lookup failed): {skipped_ctx_key}")
        print(f"Debug (TrinucModel.fit): ------------------------")
        sys.stdout.flush()

        # Calculate relative rates
        total_counts = np.sum(counts)
        if total_counts > 0 and processed_mutations > 0:
            # Normalize rates relative to the mean rate across all 192 contexts
            mean_rate = total_counts / 192.0
            # Avoid division by zero if mean_rate is somehow zero
            if mean_rate == 0:
                 print("Debug (TrinucModel.fit): ERROR - Mean rate calculated as zero. Setting rates to uniform 1.0.")
                 self.rates = np.ones((6, 16, 2))
            else:
                 self.rates = counts / mean_rate
                 print(f"Debug (TrinucModel.fit): Calculated rates. Total counts: {total_counts:.0f}")
                 print(f"Debug (TrinucModel.fit): Mean rate used for normalization: {mean_rate:.4f}")
                 print(f"Debug (TrinucModel.fit): Resulting rates - Min: {np.min(self.rates):.4f}, Max: {np.max(self.rates):.4f}, Mean: {np.mean(self.rates):.4f}")
        else:
            # If no mutations were processed, initialize with uniform rates
            print("Warning (TrinucModel.fit): No mutations matched criteria. Using uniform rates (all 1.0).")
            self.rates = np.ones((6, 16, 2))  # Default to uniform rates

        # Clear cache after fitting
        self.cache = {}
        print("Debug (TrinucModel.fit): Model fitting finished.")
        sys.stdout.flush()
        return self

    # --- get_rate method ---
    # Needs to use the same corrected context key logic
    def get_rate(self, ref: str, alt: str, context: str) -> float:
        """
        Get the mutation rate for a specific mutation (fixed context key logic)
        """
        cache_key = (ref, alt, context)
        if cache_key in self.cache:
            return self.cache[cache_key]

        # --- Input Validation ---
        if not isinstance(ref, str) or not isinstance(alt, str) or not isinstance(context, str) or \
           len(ref) != 1 or len(alt) != 1 or len(context) != 3:
            # Minimal debug print here as it's called very often
            # if debug_counters['get_rate_invalid_input'] < DEBUG_LIMIT: ...
            return 0.0

        # Handle ambiguous bases in input immediately
        if ref not in "ACGT" or alt not in "ACGT" or \
           context[0] not in "ACGTN" or context[1] not in "ACGTN" or context[2] not in "ACGTN":
            # if debug_counters['get_rate_ambiguous'] < DEBUG_LIMIT: ...
            return 0.0

        # Check context consistency (middle base should match ref)
        if context[1] != ref:
             # Attempt to fix - this might hide upstream issues if context is wrong
             # if debug_counters['get_rate_ctx_mismatch'] < DEBUG_LIMIT: ...
             context = context[0] + ref + context[2]

        # --- Normalization ---
        norm_ref, norm_alt, norm_context, strand = self._normalize_mutation(ref, alt, context)
        if strand == -1:
             return 0.0 # Normalization failed

        # --- Get Indices ---
        sub_key = f"{norm_ref}>{norm_alt}"
        if sub_key not in self.sub_to_idx:
             # if debug_counters['get_rate_sub_key'] < DEBUG_LIMIT: ...
             return 0.0  # Not a standard substitution

        sub_idx = self.sub_to_idx[sub_key]

        # Context key uses flanking bases from normalized context and 'N' in middle
        if not isinstance(norm_context, str) or len(norm_context) != 3:
             return 0.0 # Invalid normalized context

        # ***** THE FIX IS HERE *****
        ctx_key = f"{norm_context[0]}N{norm_context[2]}"
        # ***************************

        if ctx_key not in self.ctx_to_idx:
            # if debug_counters['get_rate_ctx_key'] < DEBUG_LIMIT: ...
             return 0.0  # Invalid context key

        ctx_idx = self.ctx_to_idx[ctx_key]

        # --- Get Rate ---
        rate = 0.0 # Default to 0
        # Check bounds before accessing numpy array
        if 0 <= sub_idx < self.rates.shape[0] and \
           0 <= ctx_idx < self.rates.shape[1] and \
           0 <= strand < self.rates.shape[2]:
            rate = self.rates[sub_idx, ctx_idx, strand]
        # else: # Index out of bounds (should not happen with checks)
            # if debug_counters['get_rate_index_error'] < DEBUG_LIMIT: ...

        # Store in cache
        self.cache[cache_key] = rate
        return rate
