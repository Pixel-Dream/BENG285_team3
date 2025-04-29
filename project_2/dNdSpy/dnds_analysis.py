# dnds_analysis.py - With Debugging
import os
import pandas as pd
import numpy as np
from typing import Dict, List, Tuple, Optional, Union, Any
from collections import defaultdict
from tqdm.notebook import tqdm # Use standard tqdm if not in notebook
import sys

# Import our modules (ensure these point to the debug versions)
from trinucleotide_model import TrinucleotideModel
from negative_binomial_model import NegativeBinomialModel
from dnds_calculator import DnDsCalculator
from selection_tester import SelectionTester
from mutation import Mutation, Gene, Sample, MutationDataset, ReferenceGenome, GeneAnnotation # Import necessary classes

# Define how many debug messages to print per category
DEBUG_LIMIT = 5
debug_counters = defaultdict(int)


class DnDsAnalysis:
    """
    Integrated pipeline for dN/dS analysis (with debugging prints).
    """

    def __init__(self, mutation_dataset: MutationDataset, reference_genome: ReferenceGenome,
                 gene_annotation: GeneAnnotation, covariates_df: Optional[pd.DataFrame] = None):
        """Initialize the analysis pipeline"""
        print("Debug (Analysis.__init__): Initializing DnDsAnalysis...")
        if not isinstance(mutation_dataset, MutationDataset): raise TypeError("mutation_dataset must be a MutationDataset object")
        if not isinstance(reference_genome, ReferenceGenome): raise TypeError("reference_genome must be a ReferenceGenome object")
        if not isinstance(gene_annotation, GeneAnnotation): raise TypeError("gene_annotation must be a GeneAnnotation object")

        self.dataset = mutation_dataset
        self.ref_genome = reference_genome
        self.gene_annotation = gene_annotation
        self.covariates_df = covariates_df # Optional

        print("Debug (Analysis.__init__): Initializing component models...")
        # Initialize statistical models
        self.trinuc_model = TrinucleotideModel() # Will print its own debug info
        self.nb_model = NegativeBinomialModel() # Assumes this class exists
        self.dnds_calculator = DnDsCalculator(self.trinuc_model) # Will print its own debug info
        self.selection_tester = SelectionTester() # Assumes this class exists

        # Results storage
        self.gene_results = None
        self.global_dnds = None
        print("Debug (Analysis.__init__): DnDsAnalysis initialization complete.")
        sys.stdout.flush()


    def run_analysis(self, output_dir: Optional[str] = None,
                    hypermutator_threshold: int = 500, # Note: filtering should happen before init
                    fdr_threshold: float = 0.1) -> Optional[pd.DataFrame]:
        """Run the complete dN/dS analysis pipeline (with debugging)"""
        print("\n--- Debug (Analysis.run): Starting dN/dS analysis run ---")
        sys.stdout.flush()

        # Reset debug counters for this run
        global debug_counters
        debug_counters = defaultdict(int)

        # Step 1: Use the provided dataset (filtering assumed done beforehand)
        print(f"Debug (Analysis.run): Using provided dataset with {len(self.dataset.mutations)} mutations, {len(self.dataset.samples)} samples.")
        if not self.dataset.mutations:
             print("Debug (Analysis.run): ERROR - Input mutation dataset is empty. Cannot proceed.")
             return None

        # Step 2: Fit trinucleotide model
        print("\nDebug (Analysis.run): Step 2 - Fitting trinucleotide substitution model...")
        sys.stdout.flush()
        mutation_count = len(self.dataset.mutations)
        print(f"Debug (Analysis.run): Using {mutation_count} mutations from dataset to fit model.")
        try:
            # The fit method now has detailed internal debugging
            self.trinuc_model.fit(self.dataset.mutations, self.ref_genome)
            print("Debug (Analysis.run): Trinucleotide model fitting process completed.")
            sys.stdout.flush()

            # Verify model rates after fitting
            if np.all(self.trinuc_model.rates == 1.0):
                print("Debug (Analysis.run): WARNING - Trinucleotide model rates are all 1.0 (uniform). This might indicate fitting issues or fallback.")
            elif np.all(self.trinuc_model.rates == 0.0):
                 print("Debug (Analysis.run): ERROR - Trinucleotide model rates are all 0.0! This will cause issues.")
                 # Decide how to handle this - maybe force uniform?
                 # self.trinuc_model.rates = np.ones_like(self.trinuc_model.rates)
                 # print("Debug (Analysis.run): Forcing uniform rates due to all zeros.")
            else:
                 print(f"Debug (Analysis.run): Trinuc rates seem fitted (Min={np.min(self.trinuc_model.rates):.3f}, Max={np.max(self.trinuc_model.rates):.3f}, Mean={np.mean(self.trinuc_model.rates):.3f}).")

        except Exception as e:
             print(f"Debug (Analysis.run): ERROR during trinucleotide model fitting: {e}")
             import traceback
             traceback.print_exc()
             return None # Cannot proceed without a model
        sys.stdout.flush()


        # Step 3: Process genes - Get sequences and calculate expected counts
        print("\nDebug (Analysis.run): Step 3 - Processing genes (sequences and expected counts)...")
        sys.stdout.flush()
        gene_sequences = {}
        expected_counts = {} # Stores dict: {gene_name: {'synonymous': {mut:rate,...}, 'missense': ...}}
        observed_mutations_by_gene = defaultdict(list) # Store observed mutations grouped by gene
        genes_processed = 0
        genes_skipped_no_info = 0
        genes_skipped_no_cds = 0
        genes_with_zero_expected = 0

        # Pre-group mutations by gene for faster lookup
        print("Debug (Analysis.run): Pre-grouping observed mutations by gene...")
        for mut in self.dataset.mutations:
             if mut.gene: # Ensure gene name exists
                 observed_mutations_by_gene[mut.gene].append(mut)
        print(f"Debug (Analysis.run): Grouped mutations for {len(observed_mutations_by_gene)} genes.")
        sys.stdout.flush()


        # Process each gene present in the mutation data's gene list
        # Use the genes dictionary from the MutationDataset as the primary loop source
        print(f"Debug (Analysis.run): Iterating through {len(self.dataset.genes)} genes found in the mutation data...")
        gene_iterator = tqdm(self.dataset.genes.items(), desc="Processing genes", total=len(self.dataset.genes), file=sys.stdout)

        for gene_name, gene_obj in gene_iterator:
            # A) Get gene annotation info
            gene_info = self.gene_annotation.get_gene_info(gene_name)
            if not gene_info:
                genes_skipped_no_info += 1
                if debug_counters['proc_gene_no_info'] < DEBUG_LIMIT:
                    # Use gene_iterator.write for tqdm compatibility
                    gene_iterator.write(f"  Debug (Proc.Gene): Skipping '{gene_name}' - No annotation info found.")
                    debug_counters['proc_gene_no_info'] += 1
                continue # Skip gene if no annotation

            # B) Extract coding sequence using annotation
            coding_sequence = self._get_coding_sequence(gene_name, gene_info) # This method now has debug prints
            if not coding_sequence: # Checks for length, N content happen inside _get_coding_sequence
                genes_skipped_no_cds += 1
                # Debug print happens inside _get_coding_sequence
                continue # Skip gene if CDS extraction failed

            # Store valid sequence
            gene_sequences[gene_name] = coding_sequence
            # Update coding length in the Gene object if not already set
            if gene_obj.coding_length is None:
                 gene_obj.coding_length = len(coding_sequence)


            # C) Calculate expected mutation counts using the calculator
            try:
                 # This method now has debug prints
                 gene_expected = self.dnds_calculator.calculate_expected_ns_ratio(coding_sequence)

                 # Check if expected counts are all zero for this gene
                 total_exp_gene = sum(sum(counts.values()) for counts in gene_expected.values())
                 if total_exp_gene == 0:
                      genes_with_zero_expected += 1
                      if debug_counters['proc_gene_zero_exp'] < DEBUG_LIMIT:
                           gene_iterator.write(f"  Debug (Proc.Gene): WARNING - Gene '{gene_name}' resulted in zero expected mutations. Check trinuc rates and CDS.")
                           debug_counters['proc_gene_zero_exp'] += 1
                      # Store the zero counts anyway, might be needed later
                      expected_counts[gene_name] = gene_expected
                 else:
                      expected_counts[gene_name] = gene_expected
                      # Optional: Print details for the first few successful genes
                      if debug_counters['proc_gene_success'] < DEBUG_LIMIT:
                           exp_syn = sum(gene_expected.get('synonymous', {}).values())
                           exp_mis = sum(gene_expected.get('missense', {}).values())
                           exp_non = sum(gene_expected.get('nonsense', {}).values())
                           gene_iterator.write(f"  Debug (Proc.Gene): Processed '{gene_name}'. Len={len(coding_sequence)}. ExpSyn={exp_syn:.3f}, ExpMis={exp_mis:.3f}, ExpNon={exp_non:.3f}")
                           debug_counters['proc_gene_success'] += 1

            except Exception as e:
                 if debug_counters['proc_gene_calc_error'] < DEBUG_LIMIT:
                      gene_iterator.write(f"  Debug (Proc.Gene): ERROR calculating expected counts for '{gene_name}': {e}")
                      debug_counters['proc_gene_calc_error'] += 1
                 continue # Skip gene if calculation fails

            genes_processed += 1

        print(f"\nDebug (Analysis.run): --- Gene Processing Summary ---")
        print(f"Debug (Analysis.run): Genes processed successfully: {genes_processed}")
        print(f"Debug (Analysis.run): Genes skipped (no annotation info): {genes_skipped_no_info}")
        print(f"Debug (Analysis.run): Genes skipped (CDS extraction failed): {genes_skipped_no_cds}")
        print(f"Debug (Analysis.run): Genes with zero expected mutations: {genes_with_zero_expected}")
        print(f"Debug (Analysis.run): Total expected counts stored for {len(expected_counts)} genes.")
        sys.stdout.flush()

        # Check if any expected counts were calculated
        if not expected_counts:
             print("Debug (Analysis.run): ERROR - No expected mutation counts could be calculated for any gene. Cannot proceed.")
             return None

        # Sanity check: Total expected synonymous mutations across all processed genes
        total_expected_syn_all = sum(sum(gene_expected.get('synonymous', {}).values()) for gene_expected in expected_counts.values())
        print(f"Debug (Analysis.run): Sanity Check - Total expected synonymous mutations across all genes: {total_expected_syn_all:.4f}")
        if total_expected_syn_all == 0:
            print("Debug (Analysis.run): CRITICAL WARNING - Total expected synonymous count is zero! Check trinucleotide model rates and sequence processing.")
            # Consider stopping or forcing uniform if this happens
        sys.stdout.flush()


        # Step 4: Fit negative binomial model (if covariates available)
        # Skipping this for now as covariates_df is None in the example
        if self.covariates_df is not None:
            print("\nDebug (Analysis.run): Step 4 - Fitting negative binomial model (SKIPPED - No covariates provided)...")
            # print("Fitting negative binomial model...")
            # self.nb_model.fit(self.dataset.genes, self.covariates_df) # Assumes fit method exists
        else:
            print("\nDebug (Analysis.run): Step 4 - Fitting negative binomial model (SKIPPED - No covariates provided)...")
        sys.stdout.flush()


        # Step 5: Test for selection
        print("\nDebug (Analysis.run): Step 5 - Testing for selection...")
        sys.stdout.flush()
        try:
             # Pass the pre-grouped observed mutations and the calculated expected counts
             gene_results_df = self.selection_tester.test_selection(
                 genes=self.dataset.genes, # Pass the Gene objects
                 observed_mutations=observed_mutations_by_gene, # Pass the grouped dict
                 expected_counts=expected_counts, # Pass the calculated expected counts
                 null_model='neutral', # Or other model
                 fdr_threshold=fdr_threshold
             )
             self.gene_results = gene_results_df # Store the DataFrame
             print(f"Debug (Analysis.run): Selection testing completed. Results shape: {self.gene_results.shape if self.gene_results is not None else 'None'}")
             if self.gene_results is not None and not self.gene_results.empty:
                  print(f"Debug (Analysis.run): Selection results head:\n{self.gene_results.head(DEBUG_LIMIT)}")
             elif self.gene_results is not None and self.gene_results.empty:
                  print("Debug (Analysis.run): WARNING - Selection testing returned an empty DataFrame.")
             else:
                  print("Debug (Analysis.run): WARNING - Selection testing returned None.")

        except Exception as e:
             print(f"Debug (Analysis.run): ERROR during selection testing: {e}")
             import traceback
             traceback.print_exc()
             # Decide whether to continue or stop
             self.gene_results = pd.DataFrame() # Create empty df to avoid downstream errors
        sys.stdout.flush()


        # Step 6: Calculate global dN/dS
        print("\nDebug (Analysis.run): Step 6 - Calculating global dN/dS...")
        sys.stdout.flush()
        try:
            # Use the same grouped observed mutations and calculated expected counts
            self.global_dnds = self._calculate_global_dnds(observed_mutations_by_gene, expected_counts)
            print(f"Debug (Analysis.run): Global dN/dS calculated:")
            for key, val in self.global_dnds.items():
                 print(f"  - {key}: {val:.4f}" if isinstance(val, float) else f"  - {key}: {val}")

        except Exception as e:
             print(f"Debug (Analysis.run): ERROR calculating global dN/dS: {e}")
             self.global_dnds = {} # Set empty dict
        sys.stdout.flush()


        # Step 7: Save results
        print("\nDebug (Analysis.run): Step 7 - Saving results...")
        if output_dir:
            try:
                os.makedirs(output_dir, exist_ok=True)
                print(f"Debug (Analysis.run): Output directory: {output_dir}")

                # Save gene results
                if self.gene_results is not None and not self.gene_results.empty:
                    gene_results_path = os.path.join(output_dir, 'gene_results.csv')
                    self.gene_results.to_csv(gene_results_path, index=False)
                    print(f"Debug (Analysis.run): Saved gene results to {gene_results_path}")
                else:
                    print("Debug (Analysis.run): Skipping saving gene results (None or empty).")

                # Save global dN/dS
                if self.global_dnds:
                    global_dnds_path = os.path.join(output_dir, 'global_dnds.csv')
                    # Convert dict to DataFrame before saving
                    pd.DataFrame([self.global_dnds]).to_csv(global_dnds_path, index=False)
                    print(f"Debug (Analysis.run): Saved global dN/dS to {global_dnds_path}")
                else:
                    print("Debug (Analysis.run): Skipping saving global dN/dS (empty).")

            except Exception as e:
                 print(f"Debug (Analysis.run): ERROR saving results to {output_dir}: {e}")
        else:
            print("Debug (Analysis.run): No output directory specified, skipping saving results.")
        sys.stdout.flush()

        print("\n--- Debug (Analysis.run): Analysis run method complete! ---")
        return self.gene_results


    def _get_coding_sequence(self, gene_name: str, gene_info: Dict[str, Any]) -> Optional[str]:
        """Extract and validate coding sequence for a gene (with debugging)"""
        # Check if gene info exists and contains exons
        if not gene_info or 'exons' not in gene_info or not gene_info['exons']:
            if debug_counters['cds_no_exons'] < DEBUG_LIMIT:
                 print(f"  Debug (GetCDS): Gene '{gene_name}' - No exon data found in annotation.")
                 debug_counters['cds_no_exons'] += 1
            return None

        chromosome = gene_info.get('chromosome')
        strand = gene_info.get('strand')
        exons = gene_info.get('exons', [])

        if not chromosome or not strand:
            if debug_counters['cds_no_chr_strand'] < DEBUG_LIMIT:
                 print(f"  Debug (GetCDS): Gene '{gene_name}' - Missing chromosome or strand in annotation.")
                 debug_counters['cds_no_chr_strand'] += 1
            return None

        # Sort exons by start position (important for correct assembly)
        exons.sort(key=lambda x: x[0])

        # Extract exon sequences
        exon_seqs = []
        fetch_failed = False
        for start, end in exons:
            # GTF is 1-based, fetch needs 0-based start, end is exclusive
            seq = self.ref_genome.fetch(chromosome, start - 1, end)
            if seq:
                exon_seqs.append(seq)
            else:
                # If any exon fetch fails, the CDS is incomplete/unreliable
                if debug_counters['cds_exon_fetch_fail'] < DEBUG_LIMIT:
                     print(f"  Debug (GetCDS): Gene '{gene_name}' - Failed to fetch exon sequence at {chromosome}:{start-1}-{end}.")
                     debug_counters['cds_exon_fetch_fail'] += 1
                fetch_failed = True
                break # Stop trying to fetch exons for this gene

        if fetch_failed or not exon_seqs:
             if not fetch_failed and debug_counters['cds_no_seqs_fetched'] < DEBUG_LIMIT:
                  print(f"  Debug (GetCDS): Gene '{gene_name}' - No exon sequences could be fetched.")
                  debug_counters['cds_no_seqs_fetched'] += 1
             return None # Return None if fetch failed or no sequences obtained

        # Combine exons into preliminary coding sequence
        coding_seq = ''.join(exon_seqs).upper() # Ensure uppercase

        # Reverse complement if on negative strand
        if strand == '-':
            coding_seq = self._reverse_complement(coding_seq)

        # --- Validation ---
        # 1. Check length requirement (at least one codon)
        if len(coding_seq) < 3:
            if debug_counters['cds_too_short'] < DEBUG_LIMIT:
                 print(f"  Debug (GetCDS): Gene '{gene_name}' - Assembled CDS too short (len={len(coding_seq)}). Sequence: '{coding_seq[:20]}...'")
                 debug_counters['cds_too_short'] += 1
            return None

        # 2. Ensure sequence length is a multiple of 3 (trim end if necessary)
        remainder = len(coding_seq) % 3
        if remainder > 0:
            if debug_counters['cds_trimmed'] < DEBUG_LIMIT:
                 print(f"  Debug (GetCDS): Gene '{gene_name}' - CDS length {len(coding_seq)} not multiple of 3. Trimming last {remainder} bases.")
                 debug_counters['cds_trimmed'] += 1
            coding_seq = coding_seq[:-remainder]
            # Re-check length after trimming
            if len(coding_seq) < 3:
                 if debug_counters['cds_too_short_after_trim'] < DEBUG_LIMIT:
                      print(f"  Debug (GetCDS): Gene '{gene_name}' - CDS too short after trimming (len={len(coding_seq)}).")
                      debug_counters['cds_too_short_after_trim'] += 1
                 return None


        # 3. Check for excessive ambiguous bases ('N')
        n_count = coding_seq.count('N')
        n_ratio = n_count / len(coding_seq)
        if n_ratio > 0.1: # Allow up to 10% N's
            if debug_counters['cds_too_many_n'] < DEBUG_LIMIT:
                 print(f"  Debug (GetCDS): Gene '{gene_name}' - CDS has too many 'N' bases ({n_count}/{len(coding_seq)} = {n_ratio:.2f} > 0.1).")
                 debug_counters['cds_too_many_n'] += 1
            return None

        # If all checks pass, return the validated coding sequence
        return coding_seq


    def _reverse_complement(self, sequence: str) -> str:
        """Get reverse complement of a sequence"""
        complement_map = {'A': 'T', 'C': 'G', 'G': 'C', 'T': 'A', 'N': 'N'}
        # Handle potential non-ACGTN chars gracefully by returning 'N'
        return ''.join(complement_map.get(base, 'N') for base in reversed(sequence))


    def _calculate_global_dnds(self, observed_mutations_by_gene: Dict[str, List[Mutation]],
                             expected_counts_by_gene: Dict[str, Dict[str, Dict[str, float]]]) -> Dict[str, float]:
        """Calculate global dN/dS ratios across all analyzed genes (with debugging)"""
        print("Debug (Global dN/dS): Calculating global ratios...")
        total_observed = defaultdict(int)
        total_expected = defaultdict(float)

        # Sum observed counts from the grouped dictionary
        print(f"Debug (Global dN/dS): Summing observed counts from {len(observed_mutations_by_gene)} genes...")
        genes_with_obs = 0
        for gene_name, mutations in observed_mutations_by_gene.items():
             # Only include genes for which we calculated expected counts
             if gene_name in expected_counts_by_gene:
                 genes_with_obs += 1
                 for mutation in mutations:
                     if mutation.is_synonymous(): total_observed['synonymous'] += 1
                     elif mutation.is_missense(): total_observed['missense'] += 1
                     elif mutation.is_nonsense(): total_observed['nonsense'] += 1
                     elif mutation.is_splice_site(): total_observed['splice_site'] += 1
                     # Add other categories if needed
        print(f"Debug (Global dN/dS): Summed observed counts from {genes_with_obs} relevant genes.")

        # Sum expected counts from the calculated dictionary
        print(f"Debug (Global dN/dS): Summing expected counts from {len(expected_counts_by_gene)} genes...")
        for gene_name, gene_expected in expected_counts_by_gene.items():
            total_expected['synonymous'] += sum(gene_expected.get('synonymous', {}).values())
            total_expected['missense'] += sum(gene_expected.get('missense', {}).values())
            total_expected['nonsense'] += sum(gene_expected.get('nonsense', {}).values())
            total_expected['splice_site'] += sum(gene_expected.get('splice_site', {}).values()) # Use get for safety
            # Add others if needed

        print(f"Debug (Global dN/dS): --- Totals ---")
        print(f"  Observed: Syn={total_observed['synonymous']}, Mis={total_observed['missense']}, Non={total_observed['nonsense']}, Spl={total_observed['splice_site']}")
        print(f"  Expected: Syn={total_expected['synonymous']:.3f}, Mis={total_expected['missense']:.3f}, Non={total_expected['nonsense']:.3f}, Spl={total_expected['splice_site']:.3f}")

        # Calculate dN/dS ratios safely
        dnds = {}
        n_syn_obs = total_observed['synonymous']
        e_syn_exp = total_expected['synonymous']

        if n_syn_obs > 0 and e_syn_exp > 0:
            # Calculate the synonymous mutation rate (observed/expected)
            syn_rate = n_syn_obs / e_syn_exp
            print(f"Debug (Global dN/dS): Synonymous rate (Obs/Exp): {syn_rate:.4f}")

            # Calculate non-synonymous rates and dN/dS ratios
            for mut_type in ['missense', 'nonsense', 'splice_site']:
                 n_obs = total_observed[mut_type]
                 e_exp = total_expected[mut_type]
                 if e_exp > 0:
                      non_syn_rate = n_obs / e_exp
                      dnds[mut_type] = non_syn_rate / syn_rate
                      print(f"Debug (Global dN/dS): {mut_type.capitalize()} rate: {non_syn_rate:.4f}, dN/dS_{mut_type}: {dnds[mut_type]:.4f}")
                 else:
                      dnds[mut_type] = np.nan # Assign NaN if expected is zero
                      print(f"Debug (Global dN/dS): {mut_type.capitalize()} rate: N/A (Exp=0), dN/dS_{mut_type}: NaN")


            # Calculate global non-synonymous dN/dS
            n_nonsyn_obs = total_observed['missense'] + total_observed['nonsense'] + total_observed['splice_site']
            e_nonsyn_exp = total_expected['missense'] + total_expected['nonsense'] + total_expected['splice_site']

            if e_nonsyn_exp > 0:
                global_nonsyn_rate = n_nonsyn_obs / e_nonsyn_exp
                dnds['global'] = global_nonsyn_rate / syn_rate
                print(f"Debug (Global dN/dS): Global NonSyn rate: {global_nonsyn_rate:.4f}, Global dN/dS: {dnds['global']:.4f}")
            else:
                dnds['global'] = np.nan
                print(f"Debug (Global dN/dS): Global NonSyn rate: N/A (Exp=0), Global dN/dS: NaN")

        else:
            # Cannot calculate dN/dS if observed or expected synonymous is zero
            print("Debug (Global dN/dS): WARNING - Cannot calculate dN/dS ratios (Observed Syn=0 or Expected Syn=0). Setting all to NaN.")
            dnds = {'missense': np.nan, 'nonsense': np.nan, 'splice_site': np.nan, 'global': np.nan}

        # Add raw counts for reference
        dnds.update({
            'n_syn': n_syn_obs,
            'n_mis': total_observed['missense'],
            'n_non': total_observed['nonsense'],
            'n_spl': total_observed['splice_site'],
            'exp_syn': e_syn_exp,
            'exp_mis': total_expected['missense'],
            'exp_non': total_expected['nonsense'],
            'exp_spl': total_expected['splice_site']
        })

        print("Debug (Global dN/dS): Finished calculation.")
        return dnds


    # --- Other methods (get_significant_genes, estimate_driver_mutations) ---
    # These depend on the results being generated correctly. Add debugging if needed.

    def get_significant_genes(self, fdr_threshold: float = 0.1) -> Optional[pd.DataFrame]:
        """Get genes under significant selection"""
        print(f"Debug (GetSignificant): Getting significant genes (q < {fdr_threshold})...")
        if self.gene_results is None:
            print("Debug (GetSignificant): ERROR - Analysis results not available.")
            return None
        if 'q_value' not in self.gene_results.columns:
             print("Debug (GetSignificant): ERROR - 'q_value' column missing in results.")
             return None

        significant_df = self.gene_results[self.gene_results['q_value'] < fdr_threshold].copy()
        print(f"Debug (GetSignificant): Found {len(significant_df)} significant genes.")
        return significant_df


    def estimate_driver_mutations(self, gene_list: Optional[List[str]] = None) -> Dict[str, float]:
        """Estimate the number of driver mutations (basic implementation)"""
        print("Debug (EstimateDrivers): Estimating driver mutations...")
        drivers = {'missense': 0.0, 'nonsense': 0.0, 'total': 0.0, 'drivers_per_sample': 0.0} # Initialize

        if self.gene_results is None or self.global_dnds is None:
            print("Debug (EstimateDrivers): ERROR - Analysis results or global dN/dS not available.")
            return drivers

        # Use the results DataFrame
        genes_df = self.gene_results
        if gene_list is not None:
            print(f"Debug (EstimateDrivers): Filtering results to provided list of {len(gene_list)} genes.")
            genes_df = genes_df[genes_df['gene_name'].isin(gene_list)]
            if genes_df.empty:
                 print("Debug (EstimateDrivers): WARNING - No genes from the provided list found in results.")
                 return drivers

        print(f"Debug (EstimateDrivers): Calculating excess mutations based on {len(genes_df)} genes...")

        # Sum observed and expected counts from the results DataFrame
        # Ensure columns exist and handle potential NaN values
        total_syn = genes_df['n_syn'].sum()
        total_mis = genes_df['n_mis'].sum()
        total_non = genes_df['n_non'].sum()
        # Add splice if present: total_spl = genes_df['n_spl'].sum()

        total_exp_syn = genes_df['exp_syn'].sum()
        total_exp_mis = genes_df['exp_mis'].sum()
        total_exp_non = genes_df['exp_non'].sum()
        # Add splice if present: total_exp_spl = genes_df['exp_spl'].sum()

        print(f"Debug (EstimateDrivers): Totals for selected genes - ObsSyn={total_syn}, ObsMis={total_mis}, ObsNon={total_non}")
        print(f"Debug (EstimateDrivers): Totals for selected genes - ExpSyn={total_exp_syn:.2f}, ExpMis={total_exp_mis:.2f}, ExpNon={total_exp_non:.2f}")


        # Calculate excess based on observed synonymous rate (ObsSyn / ExpSyn)
        if total_exp_syn > 0 and total_syn > 0: # Need observed syn to estimate background rate
            background_syn_rate = total_syn / total_exp_syn
            print(f"Debug (EstimateDrivers): Background synonymous rate = {background_syn_rate:.4f}")

            # Calculate expected background mutations for non-synonymous types
            expected_bg_mis = total_exp_mis * background_syn_rate
            expected_bg_non = total_exp_non * background_syn_rate
            # expected_bg_spl = total_exp_spl * background_syn_rate

            # Calculate excess (observed - expected background), floor at 0
            excess_mis = max(0.0, total_mis - expected_bg_mis)
            excess_non = max(0.0, total_non - expected_bg_non)
            # excess_spl = max(0.0, total_spl - expected_bg_spl)

            drivers['missense'] = excess_mis
            drivers['nonsense'] = excess_non
            # drivers['splice'] = excess_spl # Add if calculating

            drivers['total'] = excess_mis + excess_non # + excess_spl
            print(f"Debug (EstimateDrivers): Estimated excess - Missense={excess_mis:.2f}, Nonsense={excess_non:.2f}, Total={drivers['total']:.2f}")

        else:
             print("Debug (EstimateDrivers): WARNING - Cannot estimate background rate (Total ObsSyn or ExpSyn is zero). Driver estimate will be 0.")


        # Drivers per sample
        n_samples = len(self.dataset.samples)
        if n_samples > 0:
            drivers['drivers_per_sample'] = drivers['total'] / n_samples
            print(f"Debug (EstimateDrivers): Estimated drivers per sample = {drivers['drivers_per_sample']:.3f} ({drivers['total']:.2f} / {n_samples} samples)")
        else:
             print("Debug (EstimateDrivers): WARNING - No samples in dataset, cannot calculate drivers per sample.")
             drivers['drivers_per_sample'] = 0.0

        print("Debug (EstimateDrivers): Finished estimating drivers.")
        return drivers

