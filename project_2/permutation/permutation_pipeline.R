library(tidyverse)
library(biomaRt)
library(pbapply)
library(parallel)
source("/volume/data/BENG285/permutation_function.R")

args = commandArgs(trailingOnly=TRUE)
mutations = args[1]
#mutations = read_table("/volume/data/BENG285/TCGA.LUAD.fully_annotated.txt")
outPath = args[2]
numSims = args[3]
if(is.null(mutations)) mutations = read_table("/volume/data/BENG285/TCGA.LUAD.fully_annotated.txt")
if(is.null(outPath)) outPath = "/volume/data/BENG285/"
if(is.null(numSims)) numSims = 2000

gene_df <- mutations %>% group_by(Hugo_Symbol, Transcript_ID) %>% 
  summarise(n=n())

ensembl_mart <- useMart(biomart = "ensembl", dataset = "hsapiens_gene_ensembl")

attributes_to_fetch <- c("ensembl_transcript_id",
                         "transcript_start",
                         "transcript_end",
                         "cdna")

transcript_info <- list()
for(i in 1:ceiling(nrow(gene_df)/500)){
  transcript_info[[i]] <- tryCatch({
    message(paste("Fetching batch", i,"of",ceiling(nrow(gene_df)/500),": ",(i*500-499),"to",min(i*500,nrow(gene_df))))
    getBM(attributes = attributes_to_fetch,
          filters = "ensembl_transcript_id",
          values = gene_df$Transcript_ID[(i*500-499):min(i*500,nrow(gene_df))],
          mart = ensembl_mart)
  }, error = function(e) {
    stop(paste("An error occurred while fetching data with getBM: ", e$message))
    # Optional: Add fallback logic using archives here if needed
    # message("Trying Ensembl archive...")
    # tryCatch({ ... getBM with archive mart ...})
    return(NULL) # Return NULL if fetching fails
  })
}

n_mutation_type <- c("5'Flank","Missense_Mutation","5'UTR","3'UTR","Nonsense_Mutation",
                     "Intron","Splice_Site","Nonstop_Mutation",
                     "3'Flank","RNA","Translation_Start_Site")

trx_info_df <- do.call(rbind.data.frame, transcript_info)

gene_mutation_summ <- mutations %>% group_by(Transcript_ID,Variant_Classification) %>% summarise(n = n())

trx_info_df$n_mutation <- sapply(trx_info_df$ensembl_transcript_id,
                               \(g){
                                 sum(gene_mutation_summ$n[gene_mutation_summ$Transcript_ID == g & 
                                       gene_mutation_summ$Variant_Classification %in% n_mutation_type])
                               })
trx_info_df$s_mutation <- sapply(trx_info_df$ensembl_transcript_id,
                                 \(g){
                                   sum(gene_mutation_summ$n[gene_mutation_summ$Transcript_ID == g & 
                                                              gene_mutation_summ$Variant_Classification == "Silent"])
                                 })

trx_info_df <- trx_info_df[trx_info_df$n_mutation > 0 & trx_info_df$s_mutation > 0,]


num_cores <- 60 # Leave one core free
cl <- makeCluster(num_cores)

# Export necessary variables and load libraries on each worker node
clusterExport(cl, varlist = c("trx_info_df", "calculate_mutation_types","simulate_mutations","estimate_random_mutation_effects"), envir = environment())
clusterEvalQ(cl, {
  library(stats) # For wilcox.test
  library(dplyr) # Potentially useful within function
})


simulate_mutations_wrapper <- function(i) {
  cdna <- trx_info_df$cdna[i]
  
  tmp <- estimate_random_mutation_effects(cdna, trx_info_df$n_mutation[i]+trx_info_df$s_mutation[i], 
                                          num_rounds = numSims, dN = trx_info_df$n_mutation[i],dS = trx_info_df$s_mutation[i])
  # Return results
  return(tmp)
}

res_ls <- pblapply(1:nrow(trx_info_df), simulate_mutations_wrapper, cl = cl)

# Stop the cluster
stopCluster(cl)

permutation_df <- do.call(rbind.data.frame, res_ls)
trx_info_df$n_mutation_est <- permutation_df$average_non_synonymous
trx_info_df$s_mutation_est <- permutation_df$average_synonymous
trx_info_df <- left_join(trx_info_df,gene_df[,1:2],by = c("ensembl_transcript_id" = "Transcript_ID"))



trx_info_df$p_val <- permutation_df$quantile
trx_info_df$p_adj <- p.adjust(trx_info_df$p_val, method = "fdr")
sum(trx_info_df$p_adj<0.05)
write.csv(trx_info_df[,3:11],file.path(outPath,paste0("permutation_res_",numSims,".csv")))
write.csv(trx_info_df[trx_info_df$p_adj<0.05,3:11],file.path(outPath,paste0("permutation_res_",numSims,"_subset.csv")))