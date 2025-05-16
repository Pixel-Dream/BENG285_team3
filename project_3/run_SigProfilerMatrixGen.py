#!/usr/bin/env python3
"""
Extract de novo mutational signatures from TCGA LUAD MAF using SigProfilerExtractor.
"""
from SigProfilerMatrixGenerator import install as genInstall
from SigProfilerMatrixGenerator.scripts import SigProfilerMatrixGeneratorFunc as matGen

# 1. Define path to vcf file
def main_function():
    genInstall.install('GRCh37', rsync=False, bash=True)
    matrices = matGen.SigProfilerMatrixGeneratorFunc("project_3", "GRCh37", 
                                                    "data/vcf_filtered",
                                                    plot=True, exome=True, 
                                                    bed_file=None, chrom_based=False, 
                                                    tsb_stat=False, seqInfo=False, 
                                                    cushion=100)

if __name__=="__main__":
   main_function()
