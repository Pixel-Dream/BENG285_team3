#!/bin/bash

# Define paths
MAF_FILE="/Users/xbh0403/Desktop/25SP/BENG285/BENG285_team3/project_2/data/TCGA.LUAD.mutations.txt"
OUTPUT_FILE="/Users/xbh0403/Desktop/25SP/BENG285/BENG285_team3/project_3/output/sbs_96_matrix.all"

# Generate SBS matrix
python generate_sbs_matrix.py \
    --maf_file /Users/xbh0403/Desktop/25SP/BENG285/BENG285_team3/project_2/data/TCGA.LUAD.mutations.txt \
    --output /Users/xbh0403/Desktop/25SP/BENG285/BENG285_team3/project_3/output/sbs_96_matrix.all