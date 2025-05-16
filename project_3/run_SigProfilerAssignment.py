#!/usr/bin/env python3
"""
Assign COSMIC signatures to an SBS-96 matrix originally in CSV form.
This script reads the CSV, writes a temporary TSV, and runs cosmic_fit.
"""
import argparse
import os
import pandas as pd
from SigProfilerAssignment import Analyzer as Analyze

def main():
    parser = argparse.ArgumentParser(
        description="Convert CSV→TSV and assign COSMIC signatures with SigProfilerAssignment"
    )
    parser.add_argument(
        "-i", "--input", required=True,
        help="Path to input CSV file (comma-separated SBS96 matrix)"
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

    # 1. Read the comma-separated CSV
    df = pd.read_csv(args.input, sep=",", header=0)  # header row: MutationType + sample IDs :contentReference[oaicite:3]{index=3}

    # 2. Write out as a tab-delimited TSV
    tsv_path = os.path.splitext(args.input)[0] + ".tsv"
    df.to_csv(tsv_path, sep="\t", index=False)         # convert CSV→TSV :contentReference[oaicite:4]{index=4}

    # 3. Run SigProfilerAssignment on the TSV
    Analyze.cosmic_fit(
        samples=tsv_path,
        output=args.output,
        input_type="matrix",        # matrix of trinucleotide counts :contentReference[oaicite:5]{index=5}
        context_type=args.context_type,
        cosmic_version=args.cosmic_version,
        make_plots=True             # toggle on QC plots
    )

    print(f"Done: results in '{args.output}'")

if __name__ == "__main__":
    main()
# run wit:
# python run_SigProfilerAssignment.py \
    # -i ./W_K5.csv \
    # -o ./assignment_output
