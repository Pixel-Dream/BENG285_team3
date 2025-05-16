#!/bin/bash

# Generate SBS matrix
python matrix_generator/generate_sbs_matrix.py \
    --maf_file data/TCGA.LUAD.mutations.txt \
    --output output/sbs_96_matrix_new.tsv