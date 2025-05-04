"""
Data Preprocessing Pipeline for TCGA.LUAD.mutations.txt

This script performs the following preprocessing steps:

1. Load the original mutation annotation file (MAF format).
2. Filter to keep only mutations that passed quality control (FILTER == "PASS").
3. Remove mutations where both tumor alleles match the normal allele (likely germline).
4. Remove mutations with low tumor sequencing depth (t_depth < 10).
5. Keep only coding region mutations (specific Variant_Classification types).
6. Keep only SNPs (alleles of length 1).
7. Assign mutation weights:
    - 0.5 if one tumor allele differs from the reference allele
    - 1.0 if both tumor alleles differ from the reference allele
8. Calculate total weighted mutation counts per patient.
9. Remove hypermutator samples:
    - Samples with fewer than 500 or more than 3000 coding mutations are excluded.
10. One-hot encode the Variant_Classification column (separate mutation type counts).
11. Group data by gene (Hugo_Symbol) to summarize:
    - Total mutation weight
    - Number of each mutation type
12. Fetch gene lengths dynamically from Ensembl Biomart.
13. Merge gene length information with the gene mutation summary.
14. Save the final preprocessed data to 'result_with_gene_length.csv'.

Output:
    - A clean CSV file containing:
        * Gene symbol
        * Total weighted mutation count
        * Counts of each mutation type
        * Gene length (in base pairs)
"""

import os
import pandas as pd
from pybiomart import Server

# Step 1: Load data
df = pd.read_csv("data/TCGA.LUAD.fully_annotated.txt", sep="\t", low_memory=False)

# Step 2: Basic quality control filters
df = df[df['FILTER'] == 'PASS']

# Step 3: Remove mutations where both tumor alleles match the normal (germline)
allele_same = (
    (df['Tumor_Seq_Allele1'] == df['Match_Norm_Seq_Allele1']) &
    (df['Tumor_Seq_Allele2'] == df['Match_Norm_Seq_Allele1'])
)
df = df[~allele_same]

# Step 4: Filter by tumor depth
df = df[df['t_depth'] >= 10]

# Step 5: Keep only coding region variants
coding_variants = [
    "Missense_Mutation", "Nonsense_Mutation", "Frame_Shift_Del", "Frame_Shift_Ins",
    "Splice_Site", "Translation_Start_Site", "In_Frame_Del", "In_Frame_Ins",
    "Nonstop_Mutation"
]
df = df[df['Variant_Classification'].isin(coding_variants)]

# Step 6: Keep only SNPs (allele lengths == 1)
is_snp = (
    df['Reference_Allele'].str.len() == 1
) & (
    df['Tumor_Seq_Allele1'].str.len() == 1
) & (
    df['Tumor_Seq_Allele2'].str.len() == 1
)
df = df[is_snp]

# Step 7: Assign mutation weights
def mutation_weight(row):
    ref = row['Reference_Allele']
    allele1 = row['Tumor_Seq_Allele1']
    allele2 = row['Tumor_Seq_Allele2']
    diff1 = ref != allele1
    diff2 = ref != allele2
    return 0.5 * (diff1 + diff2)

df['mutation_weight'] = df.apply(mutation_weight, axis=1)

# Step 8: Remove hypermutator samples
mutations_per_sample = df.groupby('patient_id')['mutation_weight'].sum()
bad_samples = mutations_per_sample[(mutations_per_sample > 3000) | (mutations_per_sample < 500)].index
df = df[~df['patient_id'].isin(bad_samples)]

# Step 9: One-hot encode Variant_Classification
variant_dummies = pd.get_dummies(df['Variant_Classification'])
df = pd.concat([df, variant_dummies], axis=1)

# Step 10: Group by gene and summarize
gene_summary = df.groupby('Hugo_Symbol').agg(
    total_mutations=('mutation_weight', 'sum'),
    **{variant: (variant, 'sum') for variant in variant_dummies.columns}
).reset_index()

# Step 11: Fetch gene lengths from Ensembl Biomart
print("Fetching gene lengths from Ensembl Biomart... (please wait)")
server = Server(host='http://www.ensembl.org')
dataset = server.marts['ENSEMBL_MART_ENSEMBL'].datasets['hsapiens_gene_ensembl']

gene_info = dataset.query(attributes=['hgnc_symbol', 'start_position', 'end_position'])
gene_info['gene_length'] = abs(gene_info['Gene end (bp)'] - gene_info['Gene start (bp)']) + 1
gene_info = gene_info[['HGNC symbol', 'gene_length']]
gene_info.columns = ['Hugo_Symbol', 'gene_length']

# Step 12: Merge gene lengths into summary
final_result = gene_summary.merge(gene_info, on='Hugo_Symbol', how='left')

# Step 13: Save result
os.makedirs("results", exist_ok=True)
final_result.to_csv("results/result_with_gene_length.csv", index=False)
print("✅ Result saved to 'result_with_gene_length.csv'")
