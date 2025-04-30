import numpy as np
import os
import pandas as pd
from typing import List, Dict, Optional, Union, Any
import sys
from collections import defaultdict

# Define how many debug messages to print per category
DEBUG_LIMIT = 5
debug_counters = defaultdict(int)

class Mutation:
    """Represents a single somatic mutation event"""

    def __init__(self, patient_id: str, gene: str, chromosome: str, position: int,
                 ref_allele: str, alt_allele: str, variant_type: str,
                 variant_classification: str, **kwargs):
        self.patient_id = patient_id
        self.gene = gene
        self.chromosome = str(chromosome)
        self.position = position
        self.ref_allele = ref_allele
        self.alt_allele = alt_allele
        self.variant_type = variant_type
        self.variant_classification = variant_classification

        for key, value in kwargs.items():
            setattr(self, key, value)

    def is_synonymous(self) -> bool:
        return self.variant_classification == "Silent"

    def is_missense(self) -> bool:
        return self.variant_classification == "Missense_Mutation"

    def is_nonsense(self) -> bool:
        return self.variant_classification == "Nonsense_Mutation"

    def is_splice_site(self) -> bool:
        return self.variant_classification in ["Splice_Site", "Splice_Region"]

    def get_trinucleotide_context(self, reference_genome=None) -> Optional[str]:
        """Extract trinucleotide context for this mutation (with debugging)"""
        # Prefer context if already present in MAF
        if hasattr(self, 'context') and self.context and isinstance(self.context, str) and len(self.context) == 3:
             # Basic validation of provided context
             if all(c in "ACGTN" for c in self.context.upper()):
                # Ensure middle base matches ref allele
                if len(self.ref_allele) == 1 and self.context[1].upper() == self.ref_allele.upper():
                     return self.context.upper()
                else:
                     if debug_counters['ctx_provided_mismatch'] < DEBUG_LIMIT:
                         print(f"  Debug (Mutation.Ctx): Provided context '{self.context}' middle base mismatch ref '{self.ref_allele}' for {self}. Trying fetch.")
                         debug_counters['ctx_provided_mismatch'] += 1
             else:
                 if debug_counters['ctx_provided_invalid'] < DEBUG_LIMIT:
                     print(f"  Debug (Mutation.Ctx): Provided context '{self.context}' has invalid chars for {self}. Trying fetch.")
                     debug_counters['ctx_provided_invalid'] += 1


        # Fallback to fetching from reference genome
        if reference_genome:
            # Position in MAF is usually 1-based, fetch needs 0-based start
            # Context is ref_pos-1, ref_pos, ref_pos+1 (3 bases total)
            fetch_start = self.position - 2 # 0-based start for the base before ref
            fetch_end = self.position + 1   # 0-based end (exclusive) for the base after ref
            context = reference_genome.fetch(self.chromosome, fetch_start, fetch_end)

            if context and len(context) == 3:
                 # Basic validation
                 context = context.upper()
                 if all(c in "ACGTN" for c in context):
                     # Check if middle base matches reference allele
                     if len(self.ref_allele) == 1 and context[1] == self.ref_allele.upper():
                         return context
                     else:
                         if debug_counters['ctx_fetch_mismatch'] < DEBUG_LIMIT:
                             print(f"  Debug (Mutation.Ctx): Fetched context '{context}' middle base mismatch ref '{self.ref_allele}' for {self} at {self.chromosome}:{self.position}. Returning None.")
                             debug_counters['ctx_fetch_mismatch'] += 1
                         return None # Mismatch after fetch
                 else:
                     if debug_counters['ctx_fetch_invalid'] < DEBUG_LIMIT:
                         print(f"  Debug (Mutation.Ctx): Fetched context '{context}' has invalid chars for {self} at {self.chromosome}:{self.position}. Returning None.")
                         debug_counters['ctx_fetch_invalid'] += 1
                     return None # Invalid characters after fetch
            else:
                 # Fetch failed or returned wrong length
                 if debug_counters['ctx_fetch_fail'] < DEBUG_LIMIT:
                     print(f"  Debug (Mutation.Ctx): Fetch failed or wrong length for {self} at {self.chromosome}:{self.position} (start={fetch_start}, end={fetch_end}). Context='{context}'. Returning None.")
                     debug_counters['ctx_fetch_fail'] += 1
                 return None
        else:
             if debug_counters['ctx_no_refgenome'] < DEBUG_LIMIT:
                 print(f"  Debug (Mutation.Ctx): No reference genome provided to get context for {self}. Returning None.")
                 debug_counters['ctx_no_refgenome'] += 1
             return None # No reference genome provided

    def __str__(self) -> str:
        return f"{self.gene}:{self.chromosome}:{self.position}:{self.ref_allele}>{self.alt_allele}"

    def __repr__(self) -> str:
        return self.__str__()


class Gene:
    """Collection of mutations in a specific gene"""
    def __init__(self, gene_name: str, coding_length: Optional[int] = None,
                 gc_content: Optional[float] = None):
        self.gene_name = gene_name
        self.mutations: List[Mutation] = []
        self.coding_length = coding_length # This should be set later from annotation
        self.gc_content = gc_content
        self.expression_level = None
        self.additional_features: Dict[str, Any] = {}

    def add_mutation(self, mutation: Mutation) -> None:
        self.mutations.append(mutation)

    # Count methods remain the same...
    def count_synonymous(self) -> int: return sum(1 for m in self.mutations if m.is_synonymous())
    def count_missense(self) -> int: return sum(1 for m in self.mutations if m.is_missense())
    def count_nonsense(self) -> int: return sum(1 for m in self.mutations if m.is_nonsense())
    def count_splice_site(self) -> int: return sum(1 for m in self.mutations if m.is_splice_site())
    def count_nonsynonymous(self) -> int: return self.count_missense() + self.count_nonsense() + self.count_splice_site()

    def get_mutation_rate(self) -> Optional[float]:
        if self.coding_length and self.coding_length > 0:
            return len(self.mutations) / self.coding_length
        return None

    def __str__(self) -> str: return f"{self.gene_name}: {len(self.mutations)} mutations"
    def __repr__(self) -> str: return self.__str__()


class Sample:
    """Collection of mutations from a single patient"""
    def __init__(self, patient_id: str):
        self.patient_id = patient_id
        self.mutations: List[Mutation] = []
        self.mutation_burden = None

    def add_mutation(self, mutation: Mutation) -> None:
        self.mutations.append(mutation)

    def calculate_mutation_burden(self, exome_size: int = 30000000) -> float:
        self.mutation_burden = (len(self.mutations) / exome_size) * 1000000
        return self.mutation_burden

    def is_hypermutator(self, threshold: int = 500) -> bool:
        return len(self.mutations) > threshold

    def get_mutation_types(self) -> Dict[str, int]:
        counts = {
            'synonymous': sum(1 for m in self.mutations if m.is_synonymous()),
            'missense': sum(1 for m in self.mutations if m.is_missense()),
            'nonsense': sum(1 for m in self.mutations if m.is_nonsense()),
            'splice_site': sum(1 for m in self.mutations if m.is_splice_site())
        }
        return counts

    def __str__(self) -> str: return f"{self.patient_id}: {len(self.mutations)} mutations"
    def __repr__(self) -> str: return self.__str__()


class MutationDataset:
    """Main container for all mutation data"""
    def __init__(self):
        self.samples: Dict[str, Sample] = {}
        self.genes: Dict[str, Gene] = {}
        self.mutations: List[Mutation] = []

    def add_mutation(self, mutation: Mutation) -> None:
        """Add a mutation to the dataset and update samples and genes"""
        self.mutations.append(mutation)

        # Add to sample
        if mutation.patient_id not in self.samples:
            self.samples[mutation.patient_id] = Sample(mutation.patient_id)
        self.samples[mutation.patient_id].add_mutation(mutation)

        # Add to gene
        if mutation.gene not in self.genes:
            self.genes[mutation.gene] = Gene(mutation.gene) # Coding length added later
        self.genes[mutation.gene].add_mutation(mutation)

    def load_from_maf(self, filename: str, **filters) -> 'MutationDataset':
        """Load mutations from a MAF file (with debugging)"""
        print(f"Debug (MAF Load): Reading MAF file: {filename}")
        processed_rows = 0
        skipped_rows = 0
        skipped_key_error = 0
        skipped_value_error = 0

        try:
            # Read MAF file - consider using chunking for very large files
            maf_df = pd.read_csv(filename, sep='\t', comment='#', low_memory=False)
            print(f"Debug (MAF Load): Read {len(maf_df)} rows from MAF.")

            # --- Essential Column Check ---
            essential_cols = [
                'Hugo_Symbol', 'Chromosome', 'Start_Position',
                'Reference_Allele', 'Tumor_Seq_Allele2', # Assuming this is Alt
                'Variant_Type', 'Variant_Classification',
                'Tumor_Sample_Barcode' # Common patient ID field
            ]
            missing_cols = [col for col in essential_cols if col not in maf_df.columns]
            if missing_cols:
                 print(f"Debug (MAF Load): ERROR - MAF file is missing essential columns: {missing_cols}")
                 # Try alternative patient ID if 'Tumor_Sample_Barcode' is missing
                 if 'Tumor_Sample_Barcode' in missing_cols and 'patient_id' in maf_df.columns:
                      print("Debug (MAF Load): Found 'patient_id' column, will use that.")
                      essential_cols.remove('Tumor_Sample_Barcode')
                      essential_cols.append('patient_id')
                      missing_cols = [col for col in essential_cols if col not in maf_df.columns] # Recheck
                      if missing_cols:
                           print(f"Debug (MAF Load): ERROR - Still missing columns: {missing_cols}. Aborting load.")
                           return self # Return empty dataset
                 else:
                      print("Debug (MAF Load): Aborting MAF load.")
                      return self

            # Apply filters if specified (BEFORE iterating for efficiency)
            initial_rows = len(maf_df)
            if 'variant_classification' in filters:
                maf_df = maf_df[maf_df['Variant_Classification'].isin(filters['variant_classification'])]
                print(f"Debug (MAF Load): Filtered by Variant_Classification, {len(maf_df)} rows remain.")
            if 'variant_type' in filters:
                maf_df = maf_df[maf_df['Variant_Type'].isin(filters['variant_type'])]
                print(f"Debug (MAF Load): Filtered by Variant_Type, {len(maf_df)} rows remain.")
            # Add other filters as needed...
            if len(maf_df) < initial_rows:
                 print(f"Debug (MAF Load): Rows remaining after filtering: {len(maf_df)}")


            # Convert rows to Mutation objects and add to dataset
            print(f"Debug (MAF Load): Processing {len(maf_df)} rows...")
            for idx, row in maf_df.iterrows():
                try:
                    # Determine patient ID field
                    if 'patient_id' in maf_df.columns:
                        patient_id = row['patient_id']
                    else:
                        patient_id = row['Tumor_Sample_Barcode'] # Already checked it exists

                    # Extract required fields, handling potential NaN/None
                    gene = str(row['Hugo_Symbol']) if pd.notna(row['Hugo_Symbol']) else None
                    chrom = str(row['Chromosome']) if pd.notna(row['Chromosome']) else None
                    pos = int(row['Start_Position']) if pd.notna(row['Start_Position']) else None
                    ref = str(row['Reference_Allele']) if pd.notna(row['Reference_Allele']) else None
                    alt = str(row['Tumor_Seq_Allele2']) if pd.notna(row['Tumor_Seq_Allele2']) else None
                    v_type = str(row['Variant_Type']) if pd.notna(row['Variant_Type']) else None
                    v_class = str(row['Variant_Classification']) if pd.notna(row['Variant_Classification']) else None

                    # Basic validation of required fields
                    if not all([patient_id, gene, chrom, pos, ref, alt, v_type, v_class]):
                         raise ValueError(f"Missing essential data in row {idx}")

                    # Extract optional fields safely
                    kwargs = {}
                    optional_fields = {
                         'context': 'CONTEXT', # Common MAF field for context
                         'protein_position': 'Protein_position',
                         't_depth': 't_depth',
                         't_ref_count': 't_ref_count',
                         't_alt_count': 't_alt_count'
                    }
                    for attr_name, maf_col in optional_fields.items():
                         if maf_col in row and pd.notna(row[maf_col]):
                              # Attempt conversion for numeric fields
                              if attr_name in ['t_depth', 't_ref_count', 't_alt_count']:
                                   try:
                                        kwargs[attr_name] = int(row[maf_col])
                                   except (ValueError, TypeError):
                                        if debug_counters['maf_optional_convert_fail'] < DEBUG_LIMIT:
                                             print(f"  Debug (MAF Load): Could not convert optional field '{maf_col}' to int in row {idx}: value='{row[maf_col]}'")
                                             debug_counters['maf_optional_convert_fail'] += 1
                                        kwargs[attr_name] = None # Set to None if conversion fails
                              else:
                                   kwargs[attr_name] = str(row[maf_col]) # Store others as string

                    # Create a mutation object
                    mutation = Mutation(
                        patient_id=patient_id,
                        gene=gene,
                        chromosome=chrom,
                        position=pos,
                        ref_allele=ref,
                        alt_allele=alt,
                        variant_type=v_type,
                        variant_classification=v_class,
                        **kwargs # Pass optional fields
                    )

                    # Add to dataset
                    self.add_mutation(mutation)
                    processed_rows += 1

                except KeyError as e:
                    skipped_rows += 1
                    skipped_key_error += 1
                    if debug_counters['maf_skip_key_error'] < DEBUG_LIMIT:
                        print(f"  Debug (MAF Load): Skipping row {idx} due to KeyError: {e}")
                        debug_counters['maf_skip_key_error'] += 1
                    continue
                except (ValueError, TypeError) as e:
                    skipped_rows += 1
                    skipped_value_error += 1
                    if debug_counters['maf_skip_value_error'] < DEBUG_LIMIT:
                        print(f"  Debug (MAF Load): Skipping row {idx} due to ValueError/TypeError: {e}")
                        # print(f"    Row data: {row.to_dict()}") # Uncomment for very detailed debugging
                        debug_counters['maf_skip_value_error'] += 1
                    continue
                except Exception as e: # Catch any other unexpected errors
                     skipped_rows += 1
                     if debug_counters['maf_skip_other_error'] < DEBUG_LIMIT:
                         print(f"  Debug (MAF Load): Skipping row {idx} due to unexpected error: {e}")
                         debug_counters['maf_skip_other_error'] += 1
                     continue


            print(f"Debug (MAF Load): Finished processing. Added {processed_rows} mutations.")
            if skipped_rows > 0:
                 print(f"Debug (MAF Load): Skipped {skipped_rows} rows (KeyError: {skipped_key_error}, ValueError/TypeError: {skipped_value_error}, Other: {skipped_rows - skipped_key_error - skipped_value_error}).")

        except FileNotFoundError:
            print(f"Debug (MAF Load): ERROR - MAF file not found at {filename}")
        except Exception as e:
            print(f"Debug (MAF Load): ERROR - Unexpected error reading MAF file: {e}")
            import traceback
            traceback.print_exc()

        return self

    def filter_hypermutators(self, threshold: int = 500) -> tuple['MutationDataset', List[str]]:
        """Remove hypermutator samples from the dataset"""
        print(f"Debug (FilterHyper): Checking for hypermutators (threshold > {threshold} mutations/sample)...")
        hypermutators = []
        for sample_id, sample in self.samples.items():
             if len(sample.mutations) > threshold:
                  hypermutators.append(sample_id)

        print(f"Debug (FilterHyper): Found {len(hypermutators)} hypermutator samples: {hypermutators[:DEBUG_LIMIT]}{'...' if len(hypermutators) > DEBUG_LIMIT else ''}")

        if hypermutators:
            # Create new filtered dataset
            filtered = MutationDataset()
            kept_mutations = 0
            print(f"Debug (FilterHyper): Creating new dataset excluding hypermutators...")
            for mutation in self.mutations:
                if mutation.patient_id not in hypermutators:
                    filtered.add_mutation(mutation)
                    kept_mutations += 1
            print(f"Debug (FilterHyper): New dataset created with {kept_mutations} mutations from {len(filtered.samples)} samples.")
            return filtered, hypermutators
        else:
            print("Debug (FilterHyper): No hypermutators found. Returning original dataset.")
            return self, [] # Return self and empty list if no hypermutators

    def get_mutation_summary(self) -> Dict[str, Any]:
        """Get summary statistics for the dataset"""
        total_samples = len(self.samples)
        total_genes = len(self.genes)
        total_mutations = len(self.mutations)

        mutation_types = defaultdict(int)
        for g in self.genes.values():
             mutation_types['synonymous'] += g.count_synonymous()
             mutation_types['missense'] += g.count_missense()
             mutation_types['nonsense'] += g.count_nonsense()
             mutation_types['splice_site'] += g.count_splice_site()
             # Add counts for other types if needed

        return {
            'total_samples': total_samples,
            'total_genes': total_genes,
            'total_mutations': total_mutations,
            'mutation_types': dict(mutation_types) # Convert back to dict
        }

    def __str__(self) -> str:
        return f"MutationDataset: {len(self.mutations)} mutations in {len(self.genes)} genes from {len(self.samples)} samples"

    def __repr__(self) -> str:
        return self.__str__()


class ReferenceGenome:
    """Interface to reference genome for sequence context (with debugging)"""

    def __init__(self, fasta_file: str):
        """Initialize with a FASTA file"""
        self.fasta_path = fasta_file
        self.genome = None
        self.method = None
        self.contigs = set() # Store available contig names

        print(f"Debug (RefGenome): Initializing with FASTA: {fasta_file}")
        try:
            import pyfaidx
            print("Debug (RefGenome): Trying pyfaidx...")
            self.genome = pyfaidx.Fasta(fasta_file)
            self.method = 'pyfaidx'
            self.contigs = set(self.genome.keys())
            print(f"Debug (RefGenome): pyfaidx loaded successfully. Found {len(self.contigs)} contigs.")
            # Print a few contig names
            print(f"  Example contigs: {list(self.contigs)[:DEBUG_LIMIT]}")
        except ImportError:
            print("Debug (RefGenome): pyfaidx not found.")
            try:
                import pysam
                print("Debug (RefGenome): Trying pysam...")
                # pysam requires index (.fai file)
                if not os.path.exists(fasta_file + ".fai"):
                     print(f"Debug (RefGenome): WARNING - pysam requires an index file ({fasta_file}.fai), attempting to create...")
                     try:
                         pysam.faidx(fasta_file)
                         print(f"Debug (RefGenome): Index file created.")
                     except Exception as e:
                         print(f"Debug (RefGenome): ERROR - Failed to create index file: {e}")
                         raise ImportError("pysam requires a FASTA index (.fai), and creation failed.")

                self.genome = pysam.FastaFile(fasta_file)
                self.method = 'pysam'
                self.contigs = set(self.genome.references)
                print(f"Debug (RefGenome): pysam loaded successfully. Found {len(self.contigs)} contigs.")
                print(f"  Example contigs: {list(self.contigs)[:DEBUG_LIMIT]}")
            except ImportError:
                print("Debug (RefGenome): pysam not found either.")
                raise ImportError("ERROR: ReferenceGenome requires either pyfaidx or pysam to be installed.")
        except Exception as e:
             print(f"Debug (RefGenome): ERROR - Unexpected error during initialization: {e}")
             raise # Re-raise the exception

    def _check_chrom_format(self, chromosome: str) -> Optional[str]:
         """Check chromosome format and try alternatives"""
         # 1. Try original name
         if chromosome in self.contigs:
              return chromosome
         # 2. Try adding 'chr' if missing
         if not chromosome.startswith('chr') and f"chr{chromosome}" in self.contigs:
              if debug_counters['ref_fetch_add_chr'] < DEBUG_LIMIT:
                  print(f"  Debug (RefGenome.Fetch): Original chrom '{chromosome}' not found, using 'chr{chromosome}' instead.")
                  debug_counters['ref_fetch_add_chr'] += 1
              return f"chr{chromosome}"
         # 3. Try removing 'chr' if present
         if chromosome.startswith('chr') and chromosome[3:] in self.contigs:
              if debug_counters['ref_fetch_rem_chr'] < DEBUG_LIMIT:
                  print(f"  Debug (RefGenome.Fetch): Original chrom '{chromosome}' not found, using '{chromosome[3:]}' instead.")
                  debug_counters['ref_fetch_rem_chr'] += 1
              return chromosome[3:]
         # 4. Not found
         if debug_counters['ref_fetch_chrom_not_found'] < DEBUG_LIMIT:
             print(f"  Debug (RefGenome.Fetch): Chromosome '{chromosome}' (and variants) not found in FASTA contigs: {list(self.contigs)[:DEBUG_LIMIT]}...")
             debug_counters['ref_fetch_chrom_not_found'] += 1
         return None


    def fetch(self, chromosome: str, start: int, end: int) -> Optional[str]:
        """Fetch sequence from reference genome (with debugging)"""
        if not self.genome:
            print("Debug (RefGenome.Fetch): ERROR - Genome not loaded.")
            return None

        # Ensure start < end and start >= 0
        if start < 0 or end <= start:
             if debug_counters['ref_fetch_bad_coords'] < DEBUG_LIMIT:
                 print(f"  Debug (RefGenome.Fetch): Invalid coordinates for {chromosome}: start={start}, end={end}. Returning None.")
                 debug_counters['ref_fetch_bad_coords'] += 1
             return None


        # Check and potentially adjust chromosome name format
        valid_chromosome = self._check_chrom_format(str(chromosome)) # Ensure it's a string
        if not valid_chromosome:
             # Error printed in _check_chrom_format
             return None

        try:
            # Fetch sequence using the appropriate method
            if self.method == 'pyfaidx':
                # pyfaidx uses 0-based, half-open interval [start, end)
                sequence = str(self.genome[valid_chromosome][start:end])
                return sequence.upper() # Ensure uppercase
            else:  # pysam
                # pysam uses 0-based, half-open interval [start, end)
                sequence = self.genome.fetch(valid_chromosome, start, end)
                return sequence.upper() # Ensure uppercase

        except (KeyError, ValueError, IndexError) as e:
             # Handle cases where coordinates might be out of bounds for the chromosome
             if debug_counters['ref_fetch_coord_error'] < DEBUG_LIMIT:
                 print(f"  Debug (RefGenome.Fetch): Error fetching {valid_chromosome}:{start}-{end}. Coords likely out of bounds or other issue: {e}")
                 debug_counters['ref_fetch_coord_error'] += 1
             return None
        except Exception as e: # Catch unexpected errors
             if debug_counters['ref_fetch_unexpected'] < DEBUG_LIMIT:
                 print(f"  Debug (RefGenome.Fetch): Unexpected error fetching {valid_chromosome}:{start}-{end}: {e}")
                 debug_counters['ref_fetch_unexpected'] += 1
             return None


    def get_trinucleotide_context(self, chromosome: str, position: int) -> Optional[str]:
        """Get trinucleotide context (1-based position)"""
        # Convert 1-based position to 0-based for fetch
        # Context is pos-1, pos, pos+1
        pos0 = position - 1
        return self.fetch(chromosome, pos0 - 1, pos0 + 2) # Fetch 3 bases centered around pos0

    def get_pentanucleotide_context(self, chromosome: str, position: int) -> Optional[str]:
        """Get pentanucleotide context (1-based position)"""
        pos0 = position - 1
        return self.fetch(chromosome, pos0 - 2, pos0 + 3) # Fetch 5 bases centered around pos0

    def close(self):
        """Close the reference genome file"""
        if self.genome and hasattr(self.genome, 'close'):
            print(f"Debug (RefGenome): Closing FASTA file: {self.fasta_path}")
            self.genome.close()

# --- Gene Annotation --- (Less critical for trinuc issue, debug prints minimal)

class GeneAnnotation:
    """Interface to gene annotation for CDS information (with debugging)"""

    def __init__(self, gtf_file: Optional[str] = None):
        self.genes: Dict[str, Dict[str, Any]] = {}
        self.gtf_path = gtf_file
        if gtf_file:
            print(f"Debug (GeneAnnotation): Initializing with GTF: {gtf_file}")
            if os.path.exists(gtf_file):
                self.load_from_gtf(gtf_file)
            else:
                 print(f"Debug (GeneAnnotation): ERROR - GTF file not found: {gtf_file}")
        else:
             print("Debug (GeneAnnotation): Initialized without GTF file.")


    def _extract_attribute(self, attr_str, attr_name):
        """Extract a specific attribute from GTF attribute string"""
        for attr in attr_str.split(';'):
            attr = attr.strip()
            # Handle different quoting styles and presence/absence of space
            if attr.startswith(f'{attr_name} '):
                parts = attr.split(' ', 1) # Split only once
                if len(parts) == 2:
                    value = parts[1].strip('"')
                    return value
            elif attr.startswith(f'{attr_name}='): # Handle key=value format
                 parts = attr.split('=', 1)
                 if len(parts) == 2:
                      value = parts[1].strip('"')
                      return value
        return None

    def load_from_gtf(self, gtf_file: str) -> None:
        """Load gene annotations from GTF file"""
        print(f"Debug (GeneAnnotation): Loading gene annotations from {gtf_file}...")
        processed_genes = 0
        processed_cds_lines = 0

        try:
            # Define column names for GTF
            column_names = ['seqname', 'source', 'feature', 'start', 'end', 'score', 'strand', 'frame', 'attributes']
            # Define dtypes for faster parsing and memory efficiency
            dtypes = {
                'seqname': str, 'source': str, 'feature': str,
                'start': int, 'end': int, 'score': str, # Score often '.', handle later if needed
                'strand': str, 'frame': str, # Frame often '.', handle later if needed
                'attributes': str
            }

            # Use chunking for potentially large files
            chunk_size = 100000 # Process 100k lines at a time
            gene_data_accumulator = defaultdict(lambda: {'exons': [], 'info_set': False})

            # Read and process the GTF in chunks
            for chunk_idx, chunk in enumerate(pd.read_csv(
                    gtf_file, sep='\t', comment='#', header=None,
                    names=column_names, dtype=dtypes, chunksize=chunk_size,
                    low_memory=False # Set low_memory=False for mixed types if needed
            )):
                if (chunk_idx + 1) % 5 == 0: # Print progress every 5 chunks
                     print(f"  Debug (GeneAnnotation): Processing chunk {chunk_idx+1}...")
                     sys.stdout.flush()

                # Filter for CDS features ONLY within the chunk
                cds_chunk = chunk[chunk['feature'] == 'CDS'].copy()
                processed_cds_lines += len(cds_chunk)

                if cds_chunk.empty:
                    continue

                # Extract gene_name (or gene_id as fallback)
                # Vectorized extraction is much faster than .apply()
                cds_chunk['gene_name'] = cds_chunk['attributes'].str.extract(r'gene_name "([^"]+)"', expand=False)
                # Fallback to gene_id if gene_name is missing
                missing_name_mask = cds_chunk['gene_name'].isna()
                if missing_name_mask.any():
                     cds_chunk.loc[missing_name_mask, 'gene_name'] = cds_chunk.loc[missing_name_mask, 'attributes'].str.extract(r'gene_id "([^"]+)"', expand=False)

                # Drop rows where gene name couldn't be extracted
                cds_chunk.dropna(subset=['gene_name'], inplace=True)

                if cds_chunk.empty:
                     continue

                # Extract other info if needed (only once per gene)
                # This part is tricky with chunking, better done after aggregation
                # cds_chunk['transcript_id'] = cds_chunk['attributes'].str.extract(r'transcript_id "([^"]+)"', expand=False)

                # Aggregate exon information per gene
                for _, row in cds_chunk.iterrows():
                    gene_name = row['gene_name']
                    gene_entry = gene_data_accumulator[gene_name]

                    # Store basic info only once (from the first encountered CDS line for that gene)
                    if not gene_entry['info_set']:
                         gene_entry['chromosome'] = row['seqname']
                         gene_entry['strand'] = row['strand']
                         # Initialize coding_length
                         gene_entry['coding_length'] = 0
                         gene_entry['info_set'] = True

                    # Add exon coordinates and update coding length
                    # Ensure start <= end
                    start = min(row['start'], row['end'])
                    end = max(row['start'], row['end'])
                    gene_entry['exons'].append((start, end))
                    gene_entry['coding_length'] += (end - start + 1) # Add length of this CDS segment


            print(f"\nDebug (GeneAnnotation): Finished reading GTF. Processed {processed_cds_lines} CDS lines.")
            print(f"Debug (GeneAnnotation): Aggregating data for {len(gene_data_accumulator)} unique gene names...")

            # Finalize the gene data structure
            final_genes = {}
            for gene_name, data in gene_data_accumulator.items():
                 if data['info_set']: # Only include genes where we found CDS info
                      # Sort exons by start position
                      data['exons'].sort(key=lambda x: x[0])
                      # Remove redundant fields before storing
                      del data['info_set']
                      final_genes[gene_name] = data
                      processed_genes += 1

            self.genes = final_genes
            print(f"Debug (GeneAnnotation): Successfully processed and stored annotation for {processed_genes} genes.")

        except FileNotFoundError:
             print(f"Debug (GeneAnnotation): ERROR - GTF file not found at {gtf_file}")
        except Exception as e:
             print(f"Debug (GeneAnnotation): ERROR - Unexpected error loading GTF: {e}")
             import traceback
             traceback.print_exc()


    def get_gene_info(self, gene_name: str) -> Optional[Dict[str, Any]]:
        """Get information for a gene"""
        return self.genes.get(gene_name, None)

    def get_coding_length(self, gene_name: str) -> int:
        """Get coding length for a gene"""
        gene_info = self.get_gene_info(gene_name)
        return gene_info.get('coding_length', 0) if gene_info else 0

    def calculate_gc_content(self, gene_name: str, reference_genome) -> Optional[float]:
        """Calculate GC content for a gene's coding sequence"""
        gene_info = self.get_gene_info(gene_name)
        if not gene_info or not reference_genome:
            return None

        # Extract sequence for each exon
        sequences = []
        chromosome = gene_info.get('chromosome')
        if not chromosome: return None

        for start, end in gene_info.get('exons', []):
            # GTF is 1-based, fetch needs 0-based start
            seq = reference_genome.fetch(chromosome, start - 1, end)
            if seq:
                sequences.append(seq)
            else:
                 # If any exon fetch fails, GC content is unreliable
                 if debug_counters['gc_fetch_fail'] < DEBUG_LIMIT:
                     print(f"  Debug (GC Content): Failed to fetch exon sequence for {gene_name} at {chromosome}:{start-1}-{end}")
                     debug_counters['gc_fetch_fail'] += 1
                 return None # Cannot calculate accurately

        if not sequences:
            return None

        # Calculate GC content
        combined_seq = ''.join(sequences).upper()
        if not combined_seq: return 0.0 # Avoid division by zero

        gc_count = combined_seq.count('G') + combined_seq.count('C')
        total_len = len(combined_seq)

        return gc_count / total_len if total_len > 0 else 0.0


# --- Helper Function ---
def parse_maf_file(filename: str, **filters) -> MutationDataset:
    """Parse a MAF file and convert it to a MutationDataset"""
    print(f"Debug (ParseMAF): Creating MutationDataset and calling load_from_maf for {filename}")
    dataset = MutationDataset()
    dataset.load_from_maf(filename, **filters)
    print(f"Debug (ParseMAF): load_from_maf finished. Dataset summary: {dataset}")
    return dataset


# Example usage block removed to avoid running automatically when imported
# if __name__ == "__main__":
#     # ... example usage ...
