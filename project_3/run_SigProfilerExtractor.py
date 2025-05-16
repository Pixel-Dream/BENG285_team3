#!/usr/bin/env python3
"""
Extract de novo mutational signatures from TCGA LUAD MAF using SigProfilerExtractor.
"""
from SigProfilerExtractor import sigpro as sig

# 1. Define path to vcf file
def main_function():    
    input_data = "./data/vcf_filtered" 
    # 3. Run SigProfilerExtractor
    sig.sigProfilerExtractor(
        input_type="vcf",            # tab-delimited mutational table
        output="./SigProfilerResults/SigProfilerExtractorResults",          # output folder name
        input_data=input_data,
        reference_genome="GRCh37",      # reference genome build
        opportunity_genome="GRCh37",    # same as above for context frequencies
        exome = True,                 
        context_type="SBS96",   # contexts: SBS, DBS, ID SBS96,DBS78,ID83
        minimum_signatures=2,           # minimum candidate signatures
        maximum_signatures=20,           # maximum candidate signatures
        nmf_replicates=100,             # replicates for NMF stability
        cpu=24,                           # number of cores
    )

if __name__=="__main__":
   main_function()
