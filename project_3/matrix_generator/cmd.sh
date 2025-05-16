#!/bin/bash

# Generate SBS matrix using our code
python generate_sbs_matrix.py \
    --maf_file ../data/TCGA.LUAD.mutations.txt \
    --output ../output/sbs_96_matrix.tsv

# Generate SBS matrix using SigProfiler
python maf_to_vcf.py

python run_SigProfilerMatrixGen.py

# Run matrix comparison
python matrix_comp.py
