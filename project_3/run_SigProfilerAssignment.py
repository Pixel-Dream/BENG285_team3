#!/usr/bin/env python3
"""
Assign COSMIC signatures to an SBS-96 matrix in CSV or TSV form.
This script reads the matrix, ensures a TSV version exists, and runs cosmic_fit.
"""
import argparse
import os
import pandas as pd
from SigProfilerAssignment import Analyzer as Analyze

def main():
    parser = argparse.ArgumentParser(
        description="Convert CSV/TSV→TSV and assign COSMIC signatures with SigProfilerAssignment"
    )
    parser.add_argument(
        "-i", "--input", required=True,
        help="Path to input file (CSV or TSV SBS96 matrix)"
    )
    parser.add_argument(
        "-o", "--output", default="assignment_output",
        help="Output directory for signature assignment"
    )
    parser.add_argument(
        "--cosmic_version", type=float, default=3.3,
        help="COSMIC signature version (e.g., 3.3)"
    )
    parser.add_argument(
        "--context_type", default="96",
        help="Mutational context (default: '96')"
    )
    args = parser.parse_args()

    base, ext = os.path.splitext(args.input)
    ext = ext.lower()
    # Choose delimiter and TSV path
    if ext in (".tsv", ".txt"):
        sep_in = "\t"
        tsv_path = args.input
    else:
        sep_in = ","
        tsv_path = f"{base}.tsv"

    # Read matrix with appropriate separator
    df = pd.read_csv(args.input, sep=sep_in, header=0)

    # If original was CSV, write out TSV
    if sep_in == ",":
        df.to_csv(tsv_path, sep="\t", index=False)

    # Run SigProfilerAssignment on the TSV
    Analyze.cosmic_fit(
        samples=tsv_path,
        output=args.output,
        input_type="matrix",
        context_type=args.context_type,
        cosmic_version=args.cosmic_version,
        make_plots=True
    )

    print(f"Done: results in '{args.output}'")

if __name__ == "__main__":
    main()
# run with:
# python run_SigProfilerAssignment.py \
#     -i ./W_optim_k_4.tsv \
#     -o ./assignment_output_opt
