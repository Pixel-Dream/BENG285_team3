import os
import sys
import pandas as pd
import numpy as np
from collections import defaultdict
import argparse

def process_maf_file(maf_file):
    """
    Extract single base substitutions from a MAF file.
    
    Args:
        maf_file: Path to MAF file
        
    Returns:
        Dictionary mapping sample names to lists of mutations
    """
    # Read the MAF file
    maf_df = pd.read_csv(maf_file, sep='\t', comment='#')
    
    # Dictionary to store mutations for each sample
    mutation_lists = defaultdict(list)
    
    # Process each row in the MAF file
    for _, row in maf_df.iterrows():
        # Only process SNPs (single nucleotide variants)
        if row['Variant_Type'] == 'SNP' and row['FILTER'] == 'PASS':
            sample_name = row['patient_id']
            
            # Get the context - should be 11 bases (5 + ref + 5)
            context = row['CONTEXT']
            ref = row['Reference_Allele']
            
            # Skip if context or reference is missing
            if not context or not ref:
                continue
                
            # Ensure the reference allele is a single nucleotide
            if len(ref) != 1:
                continue
            
            # For 11-base context (standard MAF format)
            if len(context) == 11:
                # Reference should be at position 5 (0-indexed)
                center_pos = 5
                trinucleotide_context = context[center_pos-1:center_pos+2]
            else:
                # For non-standard context lengths, try to find where the reference allele is
                center_pos = len(context) // 2
                trinucleotide_context = context[center_pos-1:center_pos+2] if center_pos > 0 and center_pos < len(context) - 1 else None
            
            # Verify we have a valid trinucleotide context
            if trinucleotide_context and len(trinucleotide_context) == 3:
                alt = row['Tumor_Seq_Allele2']
                
                # Verify alt is a single nucleotide
                if len(alt) != 1:
                    continue
                
                # Store the mutation with expected middle base
                mutation_lists[sample_name].append({
                    'sample': sample_name,
                    'chrom': row['Chromosome'],
                    'pos': row['Start_Position'],
                    'ref': ref,
                    'alt': alt,
                    'context': trinucleotide_context,
                    'expected_middle': ref  # Store expected middle base for validation
                })
    
    return mutation_lists

def classify_mutation(context, ref, alt, expected_middle=None):
    """
    Classify a mutation in the SBS-96 scheme.
    
    Args:
        context: Trinucleotide context (e.g., 'ACA')
        ref: Reference base (middle of context)
        alt: Alternate base
        expected_middle: The expected middle base (to verify context)
        
    Returns:
        Mutation class in the format of "A[C>A]A"
    """
    # Convert to uppercase for consistency
    ref = ref.upper()
    alt = alt.upper()
    context = context.upper()
    
    # Verify the middle base of context matches the reference
    middle_base = context[1]
    
    # If middle base doesn't match reference, try to fix it
    if middle_base != ref:
        # If we have an expected middle base and it matches the reference,
        # something is wrong with the context
        if expected_middle and expected_middle.upper() == ref:
            # Try reverse complement in case context is on opposite strand
            complement = {'A': 'T', 'T': 'A', 'C': 'G', 'G': 'C'}
            rev_comp_context = ''.join([complement[base] for base in reversed(context)])
            
            if rev_comp_context[1] == ref:
                context = rev_comp_context
            else:
                # If still no match, can't classify this mutation
                return None
        else:
            # No expected middle or it doesn't match ref either
            return None
    
    # SBS classifications are based on pyrimidine reference bases (C, T)
    # If reference is a purine (A, G), convert to its complementary base
    purine_to_pyrimidine = {'A': 'T', 'G': 'C'}
    complement = {'A': 'T', 'T': 'A', 'C': 'G', 'G': 'C'}
    
    if ref in ['A', 'G']:
        # Need to convert to complementary strand
        new_ref = purine_to_pyrimidine[ref]
        new_alt = complement[alt]
        new_context = ''.join([complement[base] for base in reversed(context)])
        return f"{new_context[0]}[{new_ref}>{new_alt}]{new_context[2]}"
    else:
        # Already pyrimidine reference
        return f"{context[0]}[{ref}>{alt}]{context[2]}"

def generate_sbs96_mutation_types():
    """
    Generate all possible SBS-96 mutation types.
    
    Returns:
        List of all possible SBS-96 mutation types
    """
    mutation_types = []
    ref_bases = ['C', 'T']  # Pyrimidine bases
    alt_bases = ['A', 'C', 'G', 'T']
    nucleotides = ['A', 'C', 'G', 'T']
    
    for ref in ref_bases:
        for alt in alt_bases:
            if ref != alt:  # Skip identity mutations
                for upstream in nucleotides:
                    for downstream in nucleotides:
                        context = upstream + ref + downstream
                        mutation_types.append(f"{upstream}[{ref}>{alt}]{downstream}")
    
    return sorted(mutation_types)

def create_sbs_matrix(mutation_lists, matrix_type='SBS96'):
    """
    Create an SBS matrix from a list of mutations.
    
    Args:
        mutation_lists: Dictionary mapping sample names to lists of mutation dictionaries
        matrix_type: Type of SBS matrix to generate ('SBS96' is standard)
        
    Returns:
        DataFrame with mutation types as rows and samples as columns
    """
    # Generate all possible mutation types for the given matrix type
    if matrix_type == 'SBS96':
        mutation_types = generate_sbs96_mutation_types()
    else:
        raise ValueError(f"Unsupported matrix type: {matrix_type}")
    
    # Initialize counts matrix
    samples = list(mutation_lists.keys())
    counts = np.zeros((len(mutation_types), len(samples)), dtype=int)
    
    # Counters for tracking
    total_mutations = 0
    classified_mutations = 0
    
    # Fill in the matrix
    for sample_idx, sample in enumerate(samples):
        sample_mutations = mutation_lists[sample]
        for mutation in sample_mutations:
            total_mutations += 1
            context = mutation['context']
            ref = mutation['ref']
            alt = mutation['alt']
            expected_middle = mutation.get('expected_middle')
            
            try:
                mutation_class = classify_mutation(context, ref, alt, expected_middle)
                if mutation_class and mutation_class in mutation_types:
                    mutation_idx = mutation_types.index(mutation_class)
                    counts[mutation_idx, sample_idx] += 1
                    classified_mutations += 1
            except Exception as e:
                # Skip mutations that can't be classified
                continue
    
    print(f"Successfully classified {classified_mutations} out of {total_mutations} mutations")
    
    # Create DataFrame
    sbs_matrix = pd.DataFrame(counts, index=mutation_types, columns=samples)
    return sbs_matrix

def process_single_maf_file(maf_file, output_file=None):
    """
    Process a single MAF file and generate an SBS matrix.
    
    Args:
        maf_file: Path to MAF file
        output_file: Optional path to save the output matrix as CSV
        
    Returns:
        DataFrame containing the SBS matrix
    """
    print(f"Processing MAF file: {maf_file}")
    
    # Process the MAF file
    mutation_lists = process_maf_file(maf_file)
    
    # Create SBS matrix
    print("Generating SBS-96 matrix...")
    sbs_matrix = create_sbs_matrix(mutation_lists)
    
    # Save matrix if output file is specified
    if output_file:
        sbs_matrix.to_csv(output_file, sep='\t')
        print(f"SBS matrix saved to {output_file}")
    
    # Print some statistics
    print(f"Generated SBS matrix with {sbs_matrix.shape[0]} mutation types and {sbs_matrix.shape[1]} samples")
    print(f"Total mutations per sample:")
    for sample in sbs_matrix.columns:
        print(f"  {sample}: {sbs_matrix[sample].sum()} mutations")
    
    return sbs_matrix

def process_maf_directory(maf_dir, output_file=None):
    """
    Process all MAF files in a directory and generate an SBS matrix.
    
    Args:
        maf_dir: Directory containing MAF files
        output_file: Optional path to save the output matrix as CSV
        
    Returns:
        DataFrame containing the SBS matrix
    """
    # Find all MAF files
    maf_files = [f for f in os.listdir(maf_dir) if f.endswith('.maf')]
    
    if not maf_files:
        raise ValueError(f"No MAF files found in {maf_dir}")
    
    print(f"Found {len(maf_files)} MAF files to process")
    
    # Process each MAF file
    all_mutations = {}
    for maf_file in maf_files:
        print(f"Processing {maf_file}...")
        
        mutation_lists = process_maf_file(os.path.join(maf_dir, maf_file))
        all_mutations.update(mutation_lists)
    
    # Create SBS matrix
    print("Generating SBS-96 matrix...")
    sbs_matrix = create_sbs_matrix(all_mutations)
    
    # Save matrix if output file is specified
    if output_file:
        sbs_matrix.to_csv(output_file)
        print(f"SBS matrix saved to {output_file}")
    
    return sbs_matrix

def main():
    """
    Main function to parse arguments and process MAF files
    """
    parser = argparse.ArgumentParser(description='Generate an SBS-96 matrix from MAF files')
    parser.add_argument('--maf_file', help='Path to a single MAF file')
    parser.add_argument('--output', required=True, help='Output file path for the SBS matrix (CSV)')
    
    args = parser.parse_args()
    
    try:
        if args.maf_file:
            # Process a single MAF file
            if not os.path.isfile(args.maf_file):
                print(f"Error: MAF file not found: {args.maf_file}")
                return 1
            os.makedirs(os.path.dirname(args.output), exist_ok=True)
            sbs_matrix = process_single_maf_file(args.maf_file, args.output)
        else:
            print("Error: Either --maf_file or --maf_dir must be specified")
            parser.print_help()
            return 1
            
        return 0
    except Exception as e:
        print(f"Error: {e}")
        return 1

if __name__ == "__main__":
    main()