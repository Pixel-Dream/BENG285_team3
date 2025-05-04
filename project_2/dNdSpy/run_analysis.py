# main_script_debug.py
import os
import pandas as pd

# Import necessary modules
# Make sure these point to the updated files with debug prints
from data_classes import parse_maf_file, ReferenceGenome, GeneAnnotation
from dnds_analysis import DnDsAnalysis

# --- Configuration ---
MAF_FILE = "../data/TCGA.LUAD.mutations.txt"
FASTA_FILE = "../data/hg19.fa"
GTF_FILE = "../data/gencode.v19.annotation.gtf"
OUTPUT_DIR = "../results/simplified_dNdSpy_results"
HYPERMUTATOR_THRESHOLD = 500
FDR_THRESHOLD = 0.1

mutation_data = parse_maf_file(MAF_FILE)
filtered_data, hypermutators = mutation_data.filter_hypermutators(threshold=HYPERMUTATOR_THRESHOLD)
ref_genome = ReferenceGenome(FASTA_FILE)
gene_annotation = GeneAnnotation(GTF_FILE)
dnds_analysis = DnDsAnalysis(
    mutation_dataset=filtered_data,
    reference_genome=ref_genome,
    gene_annotation=gene_annotation,
    covariates_df=None
)

# --- Step 6: Run Analysis ---
print("\n--- Debug: Starting dN/dS Analysis Run ---")
try:
    # Create output directory if it doesn't exist
    os.makedirs(OUTPUT_DIR, exist_ok=True)
    print(f"Debug: Results will be saved to: {OUTPUT_DIR}")

    results = dnds_analysis.run_analysis(
        output_dir=OUTPUT_DIR,
        hypermutator_threshold=HYPERMUTATOR_THRESHOLD,
        fdr_threshold=FDR_THRESHOLD
    )

except Exception as e:
    print(f"\nDebug: ERROR during dnds_analysis.run_analysis: {e}")
    import traceback
    traceback.print_exc()
    exit()

print("\n--- Debug: Script Execution Complete ---")

# Close the reference genome file if necessary
if ref_genome and hasattr(ref_genome, 'close'):
    ref_genome.close()
    print("Debug: Closed reference genome file.")

