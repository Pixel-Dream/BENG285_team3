def count_filter_values(input_file, output_file=None):
    """
    Count filter values in a MAF file and save to JSON.
    Exactly matches the pandas query: mutations[mutations["Variant_Type"] == 'SNP']["FILTER"].value_counts()
    
    Parameters:
    -----------
    input_file : str
        Path to the input MAF file
    output_file : str, optional
        Path to save the JSON output
    """
    import pandas as pd
    import json
    from tqdm import tqdm
    
    print(f"Reading file: {input_file}")
    
    # Read the full file (or use chunksize for very large files)
    try:
        # Try with tab delimiter first (standard for MAF)
        df = pd.read_csv(input_file, sep='\t', comment='#')
    except Exception as e:
        print(f"Error with tab delimiter: {str(e)}")
        print("Trying with automatic delimiter detection...")
        df = pd.read_csv(input_file, comment='#')
    
    print(f"Total variants in file: {len(df)}")
    
    # Count unique samples
    if 'Tumor_Sample_Barcode' in df.columns:
        sample_col = 'Tumor_Sample_Barcode'
    elif 'patient_id' in df.columns:
        sample_col = 'patient_id'
    else:
        print("Warning: No patient ID column found")
        sample_col = None
        
    if sample_col:
        unique_samples = df[sample_col].nunique()
        print(f"Total unique samples: {unique_samples}")
    
    # Filter to SNP variants if Variant_Type column exists
    if 'Variant_Type' in df.columns:
        print(f"Variant types in file:")
        print(df['Variant_Type'].value_counts())
        
        # Filter to SNP only - exactly like your pandas query
        snp_df = df[df['Variant_Type'] == 'SNP']
        print(f"SNP variants: {len(snp_df)}")
    else:
        print("Warning: No Variant_Type column found")
        snp_df = df
    
    # Count filter values
    if 'FILTER' in snp_df.columns:
        filter_counts = snp_df['FILTER'].value_counts().to_dict()
        print("\nFilter values for SNP variants:")
        for filter_val, count in filter_counts.items():
            print(f"  {filter_val}: {count} ({count/len(snp_df)*100:.2f}%)")
    else:
        print("Warning: No FILTER column found")
        filter_counts = {}
    
    # Calculate non-PASS count
    pass_count = filter_counts.get('PASS', 0)
    non_pass_count = len(snp_df) - pass_count
    
    # Prepare stats output
    stats = {
        "total": {
            "samples": unique_samples if sample_col else 0,
            "variants_before_filtering": len(df),
            "snp_variants": len(snp_df),
            "pass_variants": pass_count,
            "filtered_out": {
                "non_pass": non_pass_count
            },
            "filter_values": filter_counts
        }
    }
    
    # Save to JSON if requested
    if output_file:
        with open(output_file, 'w') as f:
            json.dump(stats, f, indent=2)
        print(f"Statistics saved to {output_file}")
    
    return stats


def convert_snp_pass_to_vcfs(input_file, output_dir):
    """
    Convert SNP variants with FILTER=PASS to individual VCF files.
    Exactly matches mutations[mutations["Variant_Type"] == 'SNP']["FILTER"] == "PASS"
    
    Parameters:
    -----------
    input_file : str
        Path to the input MAF file
    output_dir : str
        Directory where individual VCF files will be saved
    """
    import pandas as pd
    import os
    import json
    from datetime import datetime
    from tqdm import tqdm
    
    # Create output directory
    os.makedirs(output_dir, exist_ok=True)
    
    print(f"Reading file: {input_file}")
    
    # Read the full file (or use chunksize for very large files)
    try:
        # Try with tab delimiter first (standard for MAF)
        df = pd.read_csv(input_file, sep='\t', comment='#')
    except Exception as e:
        print(f"Error with tab delimiter: {str(e)}")
        print("Trying with automatic delimiter detection...")
        df = pd.read_csv(input_file, comment='#')
    
    print(f"Total variants in file: {len(df)}")
    
    # Identify sample ID column
    if 'Tumor_Sample_Barcode' in df.columns:
        sample_col = 'Tumor_Sample_Barcode'
    elif 'patient_id' in df.columns:
        sample_col = 'patient_id'
    else:
        raise ValueError("No patient ID column found in file")
    
    # Get total unique samples
    unique_samples = df[sample_col].nunique()
    print(f"Total unique samples: {unique_samples}")
    
    # Apply exact pandas filters: Variant_Type == 'SNP' and FILTER == 'PASS'
    if 'Variant_Type' in df.columns:
        df = df[df['Variant_Type'] == 'SNP']
        print(f"SNP variants: {len(df)}")
    else:
        print("Warning: No Variant_Type column found, using all variants")
    
    if 'FILTER' in df.columns:
        filtered_df = df[df['FILTER'] == 'PASS']
        print(f"PASS variants: {len(filtered_df)}")
    else:
        print("Warning: No FILTER column found, using all variants")
        filtered_df = df
    
    # Count variants per sample
    sample_counts = filtered_df[sample_col].value_counts()
    print(f"Average variants per sample: {sample_counts.mean():.1f}")
    print(f"Min variants: {sample_counts.min()}, Max variants: {sample_counts.max()}")
    
    # Create stats dictionary
    filter_counts = df['FILTER'].value_counts().to_dict() if 'FILTER' in df.columns else {}
    stats = {
        "total": {
            "samples": unique_samples,
            "variants_before_filtering": len(df),
            "variants_after_filtering": len(filtered_df),
            "filtered_out": {
                "non_pass": len(df) - len(filtered_df)
            },
            "filter_values": filter_counts
        },
        "samples": {}
    }
    
    # Process each sample
    for sample_id, sample_group in tqdm(filtered_df.groupby(sample_col), desc="Samples"):
        # Skip if no variants
        if len(sample_group) == 0:
            continue
        
        # Update stats for this sample
        stats["samples"][str(sample_id)] = {
            "total_variants": len(df[df[sample_col] == sample_id]),
            "pass_variants": len(sample_group),
            "filtered_out": len(df[df[sample_col] == sample_id]) - len(sample_group)
        }
        
        # Create output filename
        safe_sample_id = str(sample_id).replace('/', '_').replace('\\', '_').replace(' ', '_')
        output_file = os.path.join(output_dir, f"{safe_sample_id}.vcf")
        
        # Create VCF file for this sample
        with open(output_file, 'w') as vcf:
            # Write VCF header
            vcf.write("##fileformat=VCFv4.2\n")
            vcf.write(f"##fileDate={datetime.now().strftime('%Y%m%d')}\n")
            vcf.write("##source=ConvertedFromMAF\n")
            vcf.write("##filter=SNP variants with FILTER=PASS only\n")
            
            # Add reference genome information if available
            if 'NCBI_Build' in sample_group.columns:
                vcf.write(f"##reference={sample_group['NCBI_Build'].iloc[0]}\n")
            
            # Add contig information
            chromosomes = sample_group['Chromosome'].unique()
            for chrom in chromosomes:
                vcf.write(f"##contig=<ID={chrom}>\n")
            
            # Write INFO fields
            for col in sample_group.columns:
                if col not in ['Chromosome', 'Start_Position', 'End_Position', sample_col, 'FILTER']:
                    vcf.write(f'##INFO=<ID={col},Number=1,Type=String,Description="{col}">\n')
            
            # Add FORMAT field for genotype
            vcf.write('##FORMAT=<ID=GT,Number=1,Type=String,Description="Genotype">\n')
            
            # Write column headers
            vcf.write(f"#CHROM\tPOS\tID\tREF\tALT\tQUAL\tFILTER\tINFO\tFORMAT\t{sample_id}\n")
            
            # Process each variant
            for _, row in sample_group.iterrows():
                # Basic fields
                chrom = row['Chromosome']
                pos = int(row['Start_Position'])
                id_field = row.get('dbSNP_RS', '.') if row.get('dbSNP_RS', '.') != '.' else '.'
                
                # REF and ALT
                ref = row['Reference_Allele']
                alt = row['Tumor_Seq_Allele2']
                
                # Handle missing values
                if pd.isna(ref) or ref == '-' or ref == '':
                    ref = 'N'
                if pd.isna(alt) or alt == '-' or alt == '':
                    alt = 'N'
                
                # Quality and filter
                qual = "."
                filter_field = "PASS"
                
                # Construct INFO field - include all other columns
                info_parts = []
                for col, value in row.items():
                    if col not in ['Chromosome', 'Start_Position', 'End_Position', sample_col, 'FILTER', 'Reference_Allele', 'Tumor_Seq_Allele2']:
                        if pd.notna(value) and value not in [".", 0, "", "-"]:
                            # Clean value for VCF format
                            value_str = str(value).replace(' ', '_').replace(';', '|').replace('=', '-')
                            info_parts.append(f"{col}={value_str}")
                
                info_field = ";".join(info_parts) if info_parts else "."
                
                # Write VCF line
                vcf.write(f"{chrom}\t{pos}\t{id_field}\t{ref}\t{alt}\t{qual}\t{filter_field}\t{info_field}\tGT\t0/1\n")
    
    # Save stats to JSON
    stats_file = os.path.join(output_dir, "filter_stats.json")
    with open(stats_file, 'w') as f:
        json.dump(stats, f, indent=2)
    
    print(f"VCF files created in {output_dir}")
    print(f"Statistics saved to {stats_file}")
    return stats


if __name__ == "__main__":    
    convert_snp_pass_to_vcfs("../data/TCGA.LUAD.mutations.txt", "../data/vcf_filtered")
