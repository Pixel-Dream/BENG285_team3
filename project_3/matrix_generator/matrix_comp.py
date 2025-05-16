import pandas as pd
import re

sig_profiler_matrix = pd.read_csv("/Users/xbh0403/Desktop/25SP/BENG285/BENG285_team3/project_3/data/vcf_filtered/output/SBS/project_3.SBS96.exome", sep="\t")
# Remove filter_status column
sig_profiler_matrix = sig_profiler_matrix.drop(columns=["filter_stats"])
# Set index to MutationType
our_matrix = pd.read_csv("/Users/xbh0403/Desktop/25SP/BENG285/BENG285_team3/project_3/output/sbs_96_matrix_our_code.tsv", sep="\t")
new_columns = our_matrix.columns.tolist()
new_columns[0] = "MutationType"
our_matrix.columns = new_columns

sig_cols = list(sig_profiler_matrix.columns)
# Extract the patient ID from the full sample ID using regex
sig_cols_short = [re.sub(r'(TCGA-\w{2}-\w{4})-.*', r'\1', col) if col != 'MutationType' else col for col in sig_cols]
sig_profiler_matrix.columns = sig_cols_short
sig_profiler_matrix.head()


import numpy as np
from scipy.stats import pearsonr, spearmanr

# First check if the columns match between the two matrices
try:
    assert sig_profiler_matrix.columns.equals(our_matrix.columns)
    print("✓ Column names match between matrices")
except AssertionError:
    print("❌ Column names do not match between matrices")
    mismatched_cols = set(sig_profiler_matrix.columns) ^ set(our_matrix.columns)
    print(f"Mismatched columns: {mismatched_cols}")

# Calculate differences and summary statistics
print("\n=== Matrix Comparison Summary ===")

# Initialize lists to store differences and counts
all_diffs = []
all_counts = []
sample_diffs = {}
col_sums = []
correlations = {}

# Compare matrices column by column (sample-wise)
for col in sig_profiler_matrix.columns[1:]:  # Skip the MutationType column
    # Extract columns from both matrices
    sig_col = sig_profiler_matrix[col]
    our_col = our_matrix[col]
    
    # Calculate differences
    diff = sig_col - our_col
    
    # Get counts from both matrices
    sig_counts = sig_col.values
    our_counts = our_col.values
    avg_counts = (sig_counts + our_counts) / 2
    
    # Store differences and counts for summary stats
    col_diffs = diff.abs().values
    all_diffs.extend(col_diffs)
    all_counts.extend(avg_counts)
    
    # Calculate sample-specific difference metric
    sample_diffs[col] = np.sum(col_diffs) / np.sum(avg_counts) if np.sum(avg_counts) > 0 else 0
    
    # Store column sums for mean calculation
    col_sums.append(np.sum(avg_counts))
    
    # Calculate correlations
    pearson_corr, _ = pearsonr(sig_counts.astype(float), our_counts.astype(float))
    spearman_corr, _ = spearmanr(sig_counts.astype(float), our_counts.astype(float))
    correlations[col] = (pearson_corr, spearman_corr)

# Calculate overall statistics
mean_diff_over_count = np.sum(all_diffs) / np.sum(all_counts) if np.sum(all_counts) > 0 else 0
median_mae = np.median(all_diffs)
max_mae = np.max(all_diffs)
mean_mutations_per_sample = np.mean(col_sums)

print(f"Mean Difference over Average Count Sum: {mean_diff_over_count:.4f}")
print(f"Overall Median Absolute Error: {median_mae:.4f}")
print(f"Maximum Absolute Error: {max_mae:.4f}")
print(f"Mean Mutations per Sample: {mean_mutations_per_sample:.2f}")

# Print samples with highest differences
print("\n=== Top 5 Samples with Highest Differences ===")
top_samples = sorted(sample_diffs.items(), key=lambda x: x[1], reverse=True)[:5]
for sample, diff_ratio in top_samples:
    print(f"{sample}: {diff_ratio:.4f}")

# Print correlation statistics
print("\n=== Correlation Statistics ===")
pearson_values = [corr[0] for corr in correlations.values()]
spearman_values = [corr[1] for corr in correlations.values()]
print(f"Mean Pearson Correlation: {np.mean(pearson_values):.4f}")
print(f"Median Pearson Correlation: {np.median(pearson_values):.4f}")
print(f"Mean Spearman Correlation: {np.mean(spearman_values):.4f}")
print(f"Median Spearman Correlation: {np.median(spearman_values):.4f}")

# Calculate and print mean column-wise correlation
print("\n=== Mean Column-wise Correlation ===")
# Calculate correlation between each mutation type across all samples
mutation_types = sig_profiler_matrix['MutationType'].values
row_correlations = []
for i in range(len(mutation_types)):
    sig_row = sig_profiler_matrix.iloc[i, 1:].values.astype(float)  # Skip MutationType column and convert to float
    our_row = our_matrix.iloc[i, 1:].values.astype(float)  # Skip MutationType column and convert to float
    row_corr, _ = pearsonr(sig_row, our_row)
    row_correlations.append(row_corr)
print(f"Mean Row-wise Pearson Correlation: {np.mean(row_correlations):.4f}")
print(f"Median Row-wise Pearson Correlation: {np.median(row_correlations):.4f}")

# Print samples with lowest correlations
print("\n=== Top 5 Samples with Lowest Correlations (Pearson) ===")
lowest_pearson = sorted(correlations.items(), key=lambda x: x[1][0])[:5]
for sample, (pearson, spearman) in lowest_pearson:
    print(f"{sample}: Pearson={pearson:.4f}, Spearman={spearman:.4f}")
