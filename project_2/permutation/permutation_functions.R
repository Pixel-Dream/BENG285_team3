#' Calculate Synonymous and Non-Synonymous Mutations in an Exon
#'
#' This function takes an exon sequence, its genomic start position, and a
#' data frame of mutations, then calculates the number of synonymous and
#' non-synonymous point mutations within that exon.
#'
#' @param exon_seq A character string representing the nucleotide sequence of the exon (e.g., "ATGCGTAGC..."). Should contain only A, T, C, G.
#' @param exon_start_pos An integer representing the genomic starting position of the exon (1-based).
#' @param mutations_df A data frame containing mutation information. It must have at least two columns:
#'   - `Position`: An integer column with the genomic position (1-based) of the mutation.
#'   - `Mutation`: A character string column describing the mutation in the format "REF>ALT" (e.g., "A>G", "T>C").
#'
#' @return A named list with two elements:
#'   - `synonymous`: The count of synonymous mutations.
#'   - `non_synonymous`: The count of non-synonymous mutations.
#'   Returns NULL and prints an error message if inputs are invalid.
#'
#' @examples
#' # Example Usage
#' exon_sequence <- "ATGCGTAGCTAGCTAGCTAGCATGC" # Example exon sequence (length 25)
#' exon_start <- 1001 # Example genomic start position
#'
#' # Example mutation data frame
#' mutations <- data.frame(
#'   Position = c(1002, 1005, 1010, 1015, 1022, 1030), # Genomic positions
#'   Mutation = c("T>C", "G>A", "T>A", "G>C", "T>G", "A>T") # REF>ALT format
#' )
#'
#' results <- calculate_mutation_types(exon_sequence, exon_start, mutations)
#' print(results)
#' # Expected output might look like (depending on the exact sequence and genetic code):
#' # $synonymous
#' # [1] 1
#' # $non_synonymous
#' # [1] 3
#' # (Mutations at 1002, 1005, 1010, 1015, 1022 are within the exon range 1001-1025)
#' # (Mutation at 1030 is outside the exon range)

calculate_mutation_types <- function(exon_seq, exon_start_pos, mutations_df) {
  
  # --- Input Validation ---
  if (!is.character(exon_seq) || length(exon_seq) != 1 || nchar(exon_seq) == 0) {
    stop("Error: 'exon_seq' must be a non-empty character string.")
  }
  if (!is.numeric(exon_start_pos) || length(exon_start_pos) != 1 || exon_start_pos <= 0 || floor(exon_start_pos) != exon_start_pos) {
    stop("Error: 'exon_start_pos' must be a single positive integer.")
  }
  if (!is.data.frame(mutations_df)) {
    stop("Error: 'mutations_df' must be a data frame.")
  }
  if (!all(c("Position", "Mutation") %in% names(mutations_df))) {
    stop("Error: 'mutations_df' must contain 'Position' and 'Mutation' columns.")
  }
  if (!is.numeric(mutations_df$Position)) {
    stop("Error: 'Position' column in 'mutations_df' must be numeric.")
  }
  if (!is.character(mutations_df$Mutation)) {
    stop("Error: 'Mutation' column in 'mutations_df' must be character.")
  }
  # Check sequence characters
  if (grepl("[^ATCG]", toupper(exon_seq))) {
    warning("Warning: 'exon_seq' contains characters other than A, T, C, G. Non-standard bases will cause errors.")
  }
  # Check mutation format (simple check)
  if (any(!grepl("^[ATCG]>[ATCG]$", toupper(mutations_df$Mutation)))) {
    warning("Warning: Some entries in 'Mutation' column do not follow the 'REF>ALT' format (e.g., 'A>G'). Malformed entries will be skipped.")
  }
  
  
  # --- Standard Genetic Code (DNA codons) ---
  genetic_code <- c(
    "TTT"="F", "TTC"="F", "TTA"="L", "TTG"="L", # Phenylalanine (F), Leucine (L)
    "TCT"="S", "TCC"="S", "TCA"="S", "TCG"="S", # Serine (S)
    "TAT"="Y", "TAC"="Y", "TAA"="*", "TAG"="*", # Tyrosine (Y), Stop (*)
    "TGT"="C", "TGC"="C", "TGA"="*", "TGG"="W", # Cysteine (C), Stop (*), Tryptophan (W)
    "CTT"="L", "CTC"="L", "CTA"="L", "CTG"="L", # Leucine (L)
    "CCT"="P", "CCC"="P", "CCA"="P", "CCG"="P", # Proline (P)
    "CAT"="H", "CAC"="H", "CAA"="Q", "CAG"="Q", # Histidine (H), Glutamine (Q)
    "CGT"="R", "CGC"="R", "CGA"="R", "CGG"="R", # Arginine (R)
    "ATT"="I", "ATC"="I", "ATA"="I", "ATG"="M", # Isoleucine (I), Methionine (M) / Start
    "ACT"="T", "ACC"="T", "ACA"="T", "ACG"="T", # Threonine (T)
    "AAT"="N", "AAC"="N", "AAA"="K", "AAG"="K", # Asparagine (N), Lysine (K)
    "AGT"="S", "AGC"="S", "AGA"="R", "AGG"="R", # Serine (S), Arginine (R)
    "GTT"="V", "GTC"="V", "GTA"="V", "GTG"="V", # Valine (V)
    "GCT"="A", "GCC"="A", "GCA"="A", "GCG"="A", # Alanine (A)
    "GAT"="D", "GAC"="D", "GAA"="E", "GAG"="E", # Aspartic Acid (D), Glutamic Acid (E)
    "GGT"="G", "GGC"="G", "GGA"="G", "GGG"="G"  # Glycine (G)
  )
  
  # --- Initialization ---
  synonymous_count <- 0
  non_synonymous_count <- 0
  exon_len <- nchar(exon_seq)
  exon_end_pos <- exon_start_pos + exon_len - 1
  exon_seq_upper <- toupper(exon_seq) # Work with uppercase
  
  # --- Process Mutations ---
  for (i in 1:nrow(mutations_df)) {
    mut_pos <- mutations_df$Position[i]
    mut_str <- toupper(mutations_df$Mutation[i])
    
    # Check if mutation is within the exon's genomic range
    if (mut_pos >= exon_start_pos && mut_pos <= exon_end_pos) {
      
      # Validate mutation format (e.g., "A>G")
      if (!grepl("^[ATCG]>[ATCG]$", mut_str)) {
        warning(paste("Skipping invalid mutation format at position", mut_pos, ":", mutations_df$Mutation[i]))
        next # Skip to next mutation
      }
      
      # Calculate position within the exon sequence (1-based)
      exon_rel_pos <- mut_pos - exon_start_pos + 1
      
      # Extract reference and alternate alleles
      ref_allele <- substr(mut_str, 1, 1)
      alt_allele <- substr(mut_str, 3, 3)
      
      # Verify the reference allele matches the exon sequence
      actual_ref <- substr(exon_seq_upper, exon_rel_pos, exon_rel_pos)
      if (actual_ref != ref_allele) {
        warning(paste("Mismatch at position", mut_pos, ": Expected ref", ref_allele, "based on mutation data, but found", actual_ref, "in exon sequence. Skipping."))
        next # Skip to next mutation
      }
      # Check if ref and alt are the same
      if (ref_allele == alt_allele) {
        # Technically a mutation wasn't specified, or it's REF>REF
        # Could classify as synonymous, or skip/warn. Let's skip for clarity.
        # warning(paste("Reference and alternate alleles are the same at position", mut_pos, ". Skipping."))
        next
      }
      
      
      # Determine the codon start position within the exon (1-based)
      # (exon_rel_pos - 1) gives 0-based position
      # floor((exon_rel_pos - 1) / 3) gives 0-based codon index
      # * 3 gives 0-based start position of the codon
      # + 1 gives 1-based start position of the codon
      codon_start_in_exon <- floor((exon_rel_pos - 1) / 3) * 3 + 1
      
      # Check if the codon is complete within the exon
      if (codon_start_in_exon + 2 > exon_len) {
        warning(paste("Mutation at position", mut_pos, "affects an incomplete codon at the end of the exon sequence. Skipping."))
        next # Skip mutation affecting incomplete codon
      }
      
      # Extract the original codon
      original_codon <- substr(exon_seq_upper, codon_start_in_exon, codon_start_in_exon + 2)
      
      # Determine the position within the codon (1, 2, or 3)
      pos_in_codon <- (exon_rel_pos - 1) %% 3 + 1
      
      # Create the mutated codon
      mutated_codon <- original_codon
      substr(mutated_codon, pos_in_codon, pos_in_codon) <- alt_allele
      
      # Translate codons to amino acids
      original_aa <- genetic_code[original_codon]
      mutated_aa <- genetic_code[mutated_codon]
      
      # Handle cases where codon is not in the genetic code table (shouldn't happen with validation)
      if (is.na(original_aa) || is.na(mutated_aa)) {
        warning(paste("Could not translate codon at position", mut_pos, "(Codons:", original_codon, "->", mutated_codon,"). Skipping."))
        next
      }
      
      # Classify mutation
      if (original_aa == mutated_aa) {
        synonymous_count <- synonymous_count + 1
      } else {
        non_synonymous_count <- non_synonymous_count + 1
      }
    } # End if mutation in exon
  } # End loop through mutations
  
  # --- Return Results ---
  return(list(synonymous = synonymous_count, non_synonymous = non_synonymous_count))
}

#' Simulate Random Point Mutations in a DNA Sequence
#'
#' Introduces a specified number of random point mutations into a DNA sequence.
#' Each mutation changes a randomly selected base to one of the other three bases.
#' Positions for mutation are chosen uniquely without replacement.
#'
#' @param dna_seq A character string representing the DNA sequence (e.g., "ATGCGTAGC...").
#'        Should ideally contain only A, T, C, G (case-insensitive).
#' @param k An integer specifying the exact number of mutations to introduce.
#'        Must be non-negative and not exceed the length of `dna_seq`.
#'
#' @return A character string representing the mutated DNA sequence.
#'         Returns the original sequence if k = 0.
#'         Stops with an error for invalid inputs.
#'
#' @examples
#' original_sequence <- "ATGCATGCATGCATGC"
#' number_of_mutations <- 3
#'
#' mutated_sequence <- simulate_mutations(original_sequence, number_of_mutations)
#' print(paste("Original:", original_sequence))
#' print(paste("Mutated: ", mutated_sequence))
#' # Example Output (will vary due to randomness):
#' # [1] "Original: ATGCATGCATGCATGC"
#' # [1] "Mutated:  ATGCATTCATGCATAC" # e.g., G->T at pos 7, G->A at pos 15, G->C at pos 11
#'
#' # Example with k=0
#' mutated_sequence_k0 <- simulate_mutations(original_sequence, 0)
#' print(paste("Mutated (k=0):", mutated_sequence_k0))
#' # [1] "Mutated (k=0): ATGCATGCATGCATGC"
#'
#' # Example with large sequence and many mutations
#' long_seq <- paste(rep("ATGC", 2500), collapse="") # Length 10000
#' system.time(mut_long_seq <- simulate_mutations(long_seq, 5000))
#' # Should be reasonably fast

simulate_mutations <- function(dna_seq, k) {
  
  # --- Input Validation ---
  if (!is.character(dna_seq) || length(dna_seq) != 1 || nchar(dna_seq) == 0) {
    stop("Error: 'dna_seq' must be a non-empty character string.")
  }
  seq_len <- nchar(dna_seq)
  if (!is.numeric(k) || length(k) != 1 || k < 0 || floor(k) != k) {
    stop("Error: 'k' must be a single non-negative integer.")
  }
  if (k > seq_len) {
    stop("Error: Number of mutations 'k' cannot exceed the sequence length.")
  }
  
  # Handle k=0 case efficiently
  if (k == 0) {
    return(dna_seq)
  }
  
  # --- Preparation ---
  # Work with uppercase and split into a character vector for modification
  dna_vec <- strsplit(toupper(dna_seq), NULL)[[1]]
  bases <- c("A", "T", "C", "G")
  
  # Check for non-standard bases
  if (any(!dna_vec %in% bases)) {
    warning("Warning: Input sequence contains characters other than A, T, C, G. Mutations will only occur at A/T/C/G positions and replace with A/T/C/G.")
    # Option: filter positions to only include standard bases if needed
    # valid_indices <- which(dna_vec %in% bases)
    # seq_len <- length(valid_indices) # update seq_len for sampling
    # if (k > seq_len) stop("Error: k exceeds number of valid A/T/C/G positions.")
    # positions_to_mutate <- sample(valid_indices, k, replace = FALSE)
  } # else { # Proceed with original indices if all bases are standard
  
  # --- Select Positions ---
  # Sample k unique positions (1-based indices) to mutate
  positions_to_mutate <- sample.int(seq_len, k, replace = FALSE)
  # } # End else for standard base check
  
  
  # --- Introduce Mutations ---
  # Loop through the selected positions and apply mutations
  for (pos in positions_to_mutate) {
    original_base <- dna_vec[pos]
    
    # Ensure the original base is one of the standard bases before mutating
    if (original_base %in% bases) {
      # Determine the possible replacement bases (excluding the original)
      replacement_options <- setdiff(bases, original_base)
      
      # Sample one replacement base randomly
      new_base <- sample(replacement_options, 1)
      
      # Update the sequence vector
      dna_vec[pos] <- new_base
    } else {
      # This part is reached if non-standard bases exist and we didn't filter indices earlier
      warning(paste("Skipping mutation at position", pos, "due to non-standard base:", original_base))
    }
  }
  
  # --- Return Result ---
  # Collapse the character vector back into a single string
  mutated_seq <- paste(dna_vec, collapse = "")
  return(mutated_seq)
}

#' Estimate Average Synonymous/Non-Synonymous Mutations from Random Mutation
#'
#' This function simulates introducing a fixed number of random mutations into
#' a given DNA sequence multiple times and estimates the average number of
#' resulting synonymous and non-synonymous changes. It relies on the
#' `simulate_mutations` and `calculate_mutation_types` functions being defined.
#'
#' @param exon_seq A character string representing the DNA sequence (e.g., coding region).
#' @param num_mutations An integer, the number of random mutations to introduce in each simulation round.
#' @param num_rounds An integer, the number of simulation rounds to perform. Default is 5000.
#'
#' @return A list containing:
#'   - `average_synonymous`: The average number of synonymous mutations per round.
#'   - `average_non_synonymous`: The average number of non-synonymous mutations per round.
#'   - `num_successful_rounds`: The number of simulation rounds that completed without errors in calculation.
#'         Returns NULL and prints an error if inputs are invalid or prerequisites are missing.
#'
#' @examples
#' # Ensure prerequisite functions are defined first:
#' # source("path/to/simulate_mutations.R")
#' # source("path/to/calculate_mutation_types.R")
#' # Or define them directly in the environment as per previous examples.
#'
#' # --- Prerequisite Functions (Example Definitions - DO NOT RUN if already loaded) ---
#' # genetic_code <- c(...) # Define genetic code as in calculate_mutation_types
#' # calculate_mutation_types <- function(...) { ... } # Definition from previous step
#' # simulate_mutations <- function(...) { ... } # Definition from previous step
#' # --- End Prerequisite Definitions ---
#'
#' # Example Usage:
#' my_exon <- "ATGCGTAGCTAGCTAGCTAGCATGC" # Length 25
#' n_mut <- 3 # Introduce 3 random mutations each time
#' n_sim <- 1000 # Perform 1000 simulation rounds (use >= 5000 for better estimate)
#'
#' # Make sure the prerequisite functions exist before running!
#' if (exists("simulate_mutations") && exists("calculate_mutation_types")) {
#'   set.seed(123) # for reproducible results
#'   estimated_effects <- estimate_random_mutation_effects(my_exon, n_mut, n_sim)
#'   print(estimated_effects)
#' } else {
#'   print("Error: Prerequisite functions 'simulate_mutations' and 'calculate_mutation_types' not found.")
#' }
#'
#' # Expected Output Structure (values depend on sequence, n_mut, n_sim, and randomness):
#' # $average_synonymous
#' # [1] 0.75 # Example value
#' # $average_non_synonymous
#' # [1] 2.25 # Example value
#' # $num_successful_rounds
#' # [1] 1000

estimate_random_mutation_effects <- function(exon_seq, num_mutations, num_rounds = 5000, dN = NULL, dS=NULL) {
  
  # --- Prerequisite Check ---
  if (!exists("simulate_mutations") || !is.function(simulate_mutations)) {
    stop("Error: Function 'simulate_mutations' is not defined or not a function.")
  }
  if (!exists("calculate_mutation_types") || !is.function(calculate_mutation_types)) {
    stop("Error: Function 'calculate_mutation_types' is not defined or not a function.")
  }
  
  # --- Input Validation ---
  if (!is.character(exon_seq) || length(exon_seq) != 1 || nchar(exon_seq) == 0) {
    stop("Error: 'exon_seq' must be a non-empty character string.")
  }
  seq_len <- nchar(exon_seq)
  if (!is.numeric(num_mutations) || length(num_mutations) != 1 || num_mutations < 0 || floor(num_mutations) != num_mutations) {
    stop("Error: 'num_mutations' must be a single non-negative integer.")
  }
  if (num_mutations > seq_len) {
    stop("Error: 'num_mutations' cannot exceed the sequence length.")
  }
  if (!is.numeric(num_rounds) || length(num_rounds) != 1 || num_rounds <= 0 || floor(num_rounds) != num_rounds) {
    stop("Error: 'num_rounds' must be a single positive integer.")
  }
  
  # --- Initialization ---
  all_syn_counts <- numeric(num_rounds)
  all_non_syn_counts <- numeric(num_rounds)
  original_vec <- strsplit(toupper(exon_seq), NULL)[[1]] # Convert once
  
  # --- Simulation Loop ---
  cat("Starting simulation for", num_rounds, "rounds...\n")
  pb <- txtProgressBar(min = 0, max = num_rounds, style = 3) # Progress bar
  
  for (i in 1:num_rounds) {
    # 1. Simulate mutations using the provided function
    # Handle potential errors within simulate_mutations if necessary, though it has its own checks
    mutated_seq <- tryCatch({
      simulate_mutations(exon_seq, num_mutations)
    }, error = function(e) {
      warning(paste("Error in simulate_mutations during round", i, ":", e$message))
      return(NULL) # Or handle differently
    })
    
    if (is.null(mutated_seq)) {
      # Mark round as failed if simulation itself failed
      all_syn_counts[i] <- NA
      all_non_syn_counts[i] <- NA
      next # Skip to the next round
    }
    
    mutated_vec <- strsplit(toupper(mutated_seq), NULL)[[1]]
    
    # 2. Identify differences (actual mutations occurred)
    diff_indices <- which(original_vec != mutated_vec)
    
    # Proceed only if there are actual differences to analyze
    if (length(diff_indices) > 0) {
      # 3. Format mutations for calculate_mutation_types
      mutation_positions <- diff_indices # Positions relative to start of exon_seq (1-based)
      ref_alleles <- original_vec[diff_indices]
      alt_alleles <- mutated_vec[diff_indices]
      mutation_strings <- paste0(ref_alleles, ">", alt_alleles)
      
      mutations_df_sim <- data.frame(
        Position = mutation_positions,
        Mutation = mutation_strings,
        stringsAsFactors = FALSE
      )
      
      # 4. Calculate mutation types using the provided function
      # Use exon_start_pos = 1 because positions are relative to exon_seq itself
      counts <- tryCatch({
        calculate_mutation_types(exon_seq, 1, mutations_df_sim)
      }, error = function(e) {
        warning(paste("Error in calculate_mutation_types during round", i, ":", e$message))
        return(NULL)
      })
      
      
      # Store results, handling potential errors from calculate_mutation_types
      if (!is.null(counts) && is.list(counts) && !is.null(counts$synonymous) && !is.null(counts$non_synonymous)) {
        all_syn_counts[i] <- counts$synonymous
        all_non_syn_counts[i] <- counts$non_synonymous
      } else {
        # Mark round as failed if calculation failed
        all_syn_counts[i] <- NA
        all_non_syn_counts[i] <- NA
        if(is.null(counts)) warning(paste("Calculation returned NULL for round", i)) # More specific warning
      }
      
    } else {
      # Case where num_mutations = 0 or (highly unlikely) mutations perfectly reverted
      all_syn_counts[i] <- 0
      all_non_syn_counts[i] <- 0
    }
    setTxtProgressBar(pb, i) # Update progress bar
  } # End simulation loop
  
  close(pb) # Close progress bar
  cat("\nSimulation finished.\n")
  
  # --- Aggregate Results ---
  successful_rounds <- sum(!is.na(all_syn_counts)) # Count non-NA rounds
  if (successful_rounds == 0 && num_rounds > 0) {
    warning("All simulation rounds failed calculation. Returning NA.")
    avg_synonymous <- NA
    avg_non_synonymous <- NA
  } else if (successful_rounds < num_rounds) {
    warning(paste("Calculation failed for", num_rounds - successful_rounds, "rounds. Averages based on successful rounds."))
    avg_synonymous <- mean(all_syn_counts, na.rm = TRUE)
    avg_non_synonymous <- mean(all_non_syn_counts, na.rm = TRUE)
  } else {
    avg_synonymous <- mean(all_syn_counts, na.rm = TRUE) # na.rm is technically not needed here but safe
    avg_non_synonymous <- mean(all_non_syn_counts, na.rm = TRUE)
  }
  if(!is.null(dN) &!is.null(dS)){
    ratio_vec <- na.omit(all_syn_counts/(all_syn_counts+all_non_syn_counts))
    ratio_ <- dS/(dN+dS)
    quantile_ <- sum(ratio_>ratio_vec)/length(ratio_vec)
  }else{
    quantile_ <- NA
  }
  
  
  # --- Return Results ---
  return(list(
    average_synonymous = avg_synonymous,
    average_non_synonymous = avg_non_synonymous,
    num_successful_rounds = successful_rounds,
    quantile = quantile_
    # Optional: Could also return sd(), or the full vectors all_syn_counts, all_non_syn_counts
  ))
}

