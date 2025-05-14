import os
import pandas as pd
import pysam
import pyfaidx
import numpy as np
from collections import defaultdict
import argparse
import sys

def extract_mutations_from_vcf(vcf_file, genome_ref_fasta, sample_name=None):
    """
    Extract single base substitutions from a VCF file using pysam.
    
    Args:
        vcf_file: Path to VCF file
        genome_ref_fasta: Path to reference genome fasta file
        sample_name: Optional name for the sample (defaults to filename)
        
    Returns:
        List of mutations with context information
    """
    # Load reference genome
    genome = pyfaidx.Fasta(genome_ref_fasta)
    
    # If sample name not provided, use filename without extension
    if sample_name is None:
        sample_name = os.path.splitext(os.path.basename(vcf_file))[0]
    
    # Parse VCF file with pysam
    vcf_reader = pysam.VariantFile(vcf_file, 'r')
    mutations = []
    
    for record in vcf_reader:
        # Get basics: CHROM, POS, REF, ALT
        chrom = record.chrom
        pos = record.pos  # 1-based position in pysam
        ref = record.ref
        
        # Only process SNVs (single nucleotide variants)
        if len(ref) == 1 and len(record.alts) > 0 and len(record.alts[0]) == 1:
            alt = record.alts[0]
            
            # Get trinucleotide context (the nucleotide before and after the mutation)
            try:
                # Get the +/- 1 nucleotide context (adjust for 0-based indexing)
                context = str(genome[chrom][pos-2:pos+1].seq).upper()
                if len(context) == 3:
                    mutations.append({
                        'sample': sample_name,
                        'chrom': chrom,
                        'pos': pos,
                        'ref': ref,
                        'alt': alt,
                        'context': context
                    })
            except:
                # Skip if we can't get proper context (e.g., at chromosome boundaries)
                continue
    
    return mutations

def classify_mutation(context, ref, alt):
    """
    Classify a mutation in the SBS-96 scheme.
    
    Args:
        context: Trinucleotide context (e.g., 'ACA')
        ref: Reference base (middle of context)
        alt: Alternate base
        
    Returns:
        Mutation class in the format of "C>A:ACA"
    """
    # Verify the middle base of context matches the reference
    assert context[1] == ref, f"Reference base mismatch: {context[1]} vs {ref}"
    
    # SBS classifications are based on pyrimidine reference bases (C, T)
    # If reference is a purine (A, G), convert to its complementary base
    purine_to_pyrimidine = {'A': 'T', 'G': 'C'}
    complement = {'A': 'T', 'T': 'A', 'C': 'G', 'G': 'C'}
    
    if ref in ['A', 'G']:
        # Need to convert to complementary strand
        new_ref = purine_to_pyrimidine[ref]
        new_alt = complement[alt]
        new_context = ''.join([complement[base] for base in reversed(context)])
        return f"{new_ref}>{new_alt}:{new_context}"
    else:
        # Already pyrimidine reference
        return f"{ref}>{alt}:{context}"

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
                        mutation_types.append(f"{ref}>{alt}:{context}")
    
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
    
    # Fill in the matrix
    for sample_idx, sample in enumerate(samples):
        sample_mutations = mutation_lists[sample]
        for mutation in sample_mutations:
            context = mutation['context']
            ref = mutation['ref']
            alt = mutation['alt']
            
            mutation_class = classify_mutation(context, ref, alt)
            if mutation_class in mutation_types:
                mutation_idx = mutation_types.index(mutation_class)
                counts[mutation_idx, sample_idx] += 1
    
    # Create DataFrame
    sbs_matrix = pd.DataFrame(counts, index=mutation_types, columns=samples)
    return sbs_matrix

def process_vcf_directory(vcf_dir, genome_ref_fasta, output_file=None):
    """
    Process all VCF files in a directory and generate an SBS matrix.
    
    Args:
        vcf_dir: Directory containing VCF files
        genome_ref_fasta: Path to reference genome fasta file
        output_file: Optional path to save the output matrix as CSV
        
    Returns:
        DataFrame containing the SBS matrix
    """
    # Find all VCF files
    vcf_files = [f for f in os.listdir(vcf_dir) if f.endswith('.vcf') or f.endswith('.vcf.gz')]
    
    if not vcf_files:
        raise ValueError(f"No VCF files found in {vcf_dir}")
    
    print(f"Found {len(vcf_files)} VCF files to process")
    
    # Process each VCF file
    mutation_lists = {}
    for vcf_file in vcf_files:
        sample_name = os.path.splitext(vcf_file)[0].replace('.vcf', '')  # Handle both .vcf and .vcf.gz
        print(f"Processing {vcf_file} as sample {sample_name}...")
        
        mutations = extract_mutations_from_vcf(
            os.path.join(vcf_dir, vcf_file),
            genome_ref_fasta,
            sample_name
        )
        
        mutation_lists[sample_name] = mutations
        print(f"  Found {len(mutations)} SNVs in {sample_name}")
    
    # Create SBS matrix
    print("Generating SBS-96 matrix...")
    sbs_matrix = create_sbs_matrix(mutation_lists)
    
    # Save matrix if output file is specified
    if output_file:
        sbs_matrix.to_csv(output_file)
        print(f"SBS matrix saved to {output_file}")
    
    return sbs_matrix


def main():
    parser = argparse.ArgumentParser(description='Generate an SBS-96 matrix from VCF files')
    parser.add_argument('--vcf_dir', required=True, help='Directory containing VCF files')
    parser.add_argument('--ref_genome', required=True, help='Path to reference genome FASTA file')
    parser.add_argument('--output', required=True, help='Output file path for the SBS matrix (CSV)')
    
    args = parser.parse_args()
    
    # Process all VCF files and generate the SBS matrix
    try:
        sbs_matrix = process_vcf_directory(args.vcf_dir, args.ref_genome, args.output)
        
        # Print some statistics
        print(f"Generated SBS matrix with {sbs_matrix.shape[0]} mutation types and {sbs_matrix.shape[1]} samples")
        print(f"Total mutations per sample:")
        for sample in sbs_matrix.columns:
            print(f"  {sample}: {sbs_matrix[sample].sum()} mutations")
            
    except Exception as e:
        print(f"Error: {e}", file=sys.stderr)
        return 1
    
    return 0

if __name__ == "__main__":
    sys.exit(main())