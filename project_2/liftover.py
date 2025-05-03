#!/usr/bin/env python3

# Run with: python liftover.py -i data/TCGA.LUAD.mutations.txt -c data/hg19ToHg38.over.chain.gz -r data/GRCh38_full_analysis_set_plus_decoy_hla.fa -o data/TCGA.LUAD.mutations.hg38.txt
import argparse
import pandas as pd
from pyliftover import LiftOver
import pysam

def main():
    parser = argparse.ArgumentParser(
        description="LiftOver TCGA MAF from hg19 to hg38 with allele fetching"
    )
    parser.add_argument(
        "-i", "--input", required=True,
        help="Input MAF file (GRCh37 coordinates)"
    )
    parser.add_argument(
        "-c", "--chain", required=True,
        help="UCSC chain file (e.g., hg19ToHg38.over.chain.gz)"
    )
    parser.add_argument(
        "-r", "--reference", required=True,
        help="Indexed hg38 FASTA file"
    )
    parser.add_argument(
        "-o", "--output", required=True,
        help="Output MAF file with hg38 columns"
    )
    args = parser.parse_args()

    # Load the original MAF (chromosome as string) :contentReference[oaicite:5]{index=5}
    maf = pd.read_csv(
        args.input, sep="\t", comment="#", dtype={"Chromosome": str},
        low_memory=False
    )

    # Initialize liftover and FASTA reader :contentReference[oaicite:6]{index=6}
    lo = LiftOver(args.chain)
    fasta = pysam.FastaFile(args.reference)

    def liftover_variant(row):
        # Ensure 'chr'-prefixed chromosome for LiftOver :contentReference[oaicite:7]{index=7}
        chrom = row["Chromosome"]
        if not chrom.startswith("chr"):
            chrom = "chr" + chrom
        start = int(row["Start_Position"])
        end = int(row["End_Position"])
        # Convert 0-based start coordinate :contentReference[oaicite:8]{index=8}
        conv = lo.convert_coordinate(chrom, start - 1)
        if not conv:
            # Unmapped variants yield NaNs
            return pd.Series({
                "hg38_Chromosome": None,
                "hg38_Start_Position": None,
                "hg38_End_Position": None,
                "hg38_Reference_Allele": None,
                "hg38_Tumor_Seq_Allele1": None,
                "hg38_Tumor_Seq_Allele2": None
            })
        new_chrom, new_pos0, strand, _ = conv[0]
        new_start = new_pos0 + 1
        length = end - start + 1
        new_end = new_start + length - 1
        # Fetch the reference allele from hg38 FASTA :contentReference[oaicite:9]{index=9}
        try:
            ref_allele = fasta.fetch(new_chrom, new_start - 1, new_end)
        except:
            ref_allele = None
        return pd.Series({
            "hg38_Chromosome": new_chrom.replace("chr", ""),
            "hg38_Start_Position": new_start,
            "hg38_End_Position": new_end,
            "hg38_Reference_Allele": ref_allele,
            "hg38_Tumor_Seq_Allele1": row["Tumor_Seq_Allele1"],
            "hg38_Tumor_Seq_Allele2": row["Tumor_Seq_Allele2"]
        })

    # Apply liftover row-wise :contentReference[oaicite:10]{index=10}  
    hg38_cols = maf.apply(liftover_variant, axis=1)
    # Concatenate original and new columns :contentReference[oaicite:11]{index=11}  
    maf_hg38 = pd.concat([maf, hg38_cols], axis=1)
    # Write output with tab delimiter
    maf_hg38.to_csv(args.output, sep="\t", index=False)
    print(f"Output written to {args.output}")

if __name__ == "__main__":
    main()
