setwd("/Users/zhouhaowen/Desktop/UCSD/BENG_285")
library(tidyverse)
# Load data
mutation_df <- read_tsv("TCGA.LUAD.mutations.txt")

data_mat <- read_tsv("TCGA.LUAD.expression.txt")

expr_mat <- as.matrix(data_mat[,-c(1:2)])
meta_df <- data.frame(patient_id = data_mat$patient_id,
                      #cases.submitter_id = data_mat$patient_id,
                      sample_id = data_mat$sample_id)

row.names(expr_mat) <- as.character(meta_df$patient_id) 
colnames(expr_mat) <- sapply(colnames(expr_mat),
                             \(g){
                               tmp <- str_split(g,"\\|") %>% unlist()
                               if(tmp[1] == "?") paste0(tmp[1],tmp[2])
                               else tmp[1]
                             })
expr_mat <- expr_mat[,!str_detect(colnames(expr_mat),"\\?")]

# gene filtration
gene_df <- data.frame(gene = colnames(expr_mat),
                      expr_sum = colSums(expr_mat))
# reads < 10
expr_mat <- expr_mat[,gene_df$expr_sum >= 50]
row.names(expr_mat) <- meta_df$sample_id
# sample filtration
meta_df[["total_expr"]] <- rowSums(expr_mat)

# using course dataset
LUAD_meta <- read_tsv("TCGA.LUAD.metadata.txt")
merged_df <- left_join(meta_df, LUAD_meta, by = "patient_id", 
                       relationship = "many-to-one", multiple = "first")
row.names(merged_df) <- merged_df$sample_id


# export normalized data
require(sceasy)
require(Seurat)

obj <- CreateSeuratObject(counts = t(expr_mat), meta.data = merged_df)
obj <- NormalizeData(obj) %>% FindVariableFeatures(selection.method = "vst", nfeatures = 2000) %>% 
  ScaleData(features = rownames(.))

# find markers
Idents(obj) <- "tumor_status"
markers <- FindMarkers(obj, ident.1 = "WITH TUMOR", ident.2 = "TUMOR FREE",logfc.threshold = 0.0)
library(EnhancedVolcano)

volcano_obj <- data.frame(
  "gene" = row.names(markers),
  "logFC" = markers$avg_log2FC,
  "P.Val" = markers$p_val
)

EnhancedVolcano(volcano_obj,
                lab = volcano_obj$gene,
                x = 'logFC',
                pCutoff = 0.05,
                y = 'P.Val',
                drawConnectors = TRUE,
                legendPosition = 'right',
                labSize = 3.0,
                FCcutoff = 1,
                legendLabSize = 12,
                title = "Volcano Plot of LUAD",
                subtitle = "p-value cutoff = 0.05, log fold change cutoff=1",
                ylim = c(0,3)
)

#######
# use TCGA raw data

library(TCGAbiolinks)

query <- GDCquery(
  project      = "TCGA-LUAD",
  data.category= "Transcriptome Profiling",
  data.type    = "Gene Expression Quantification",
  workflow.type= "STAR - Counts",      # HTSeq pipeline was retired in DR‑42
  sample.type  = c("Primary Tumor","Solid Tissue Normal") # optional
)

GDCdownload(query, method = "api")      # or method="client" for the DTT
se <- GDCprepare(query, summarizedExperiment = TRUE)

se_meta <- as.data.frame(se@colData)
se_expr <- se@assays@data@listData[["tpm_unstrand"]]

colnames(se_expr) <- se_meta$barcode
row.names(se_expr) <- se@rowRanges@elementMetadata@listData[["gene_name"]]

obj <- CreateSeuratObject(counts = se_expr, meta.data = se_meta)
obj <- NormalizeData(obj) %>% FindVariableFeatures(selection.method = "vst", nfeatures = 2000) %>% 
  ScaleData(features = rownames(.))

obj <- RunPCA(obj,npcs = 50)
ElbowPlot(obj)

obj <- RunUMAP(obj,dims = 1:6)

DimPlot(obj,reduction = "pca", group.by = "tissue_type")
DimPlot(obj,reduction = "umap", group.by = "tissue_type")
DimPlot(obj,reduction = "umap", group.by = "gender")

sceasy::convertFormat(obj, from = "seurat", to = "anndata", assay = "RNA", 
                      main_layer = "scale.data", transfer_layers="data", outFile = "LUAD_TPM_normalized.h5ad")

Idents(obj) <- "tissue_type"
markers <- FindMarkers(obj, ident.1 = "Tumor", ident.2 = "Normal",logfc.threshold = 0.1)

deseq_expr <- round(t(expr_mat)*10)
deseq_meta <- merged_df
deseq_meta$tumor_status[is.na(deseq_meta$tumor_status)] <- "NA"
deseq_expr <- deseq_expr[,!deseq_meta$tumor_status %in% c("[Discrepancy]","NA")]
deseq_meta <- deseq_meta[!deseq_meta$tumor_status %in% c("[Discrepancy]","NA"),]

dds <- DESeqDataSetFromMatrix(countData = se@assays@data@listData[["unstranded"]],
                              colData = se_meta,
                              design= ~ gender + tissue_type)

dds <- DESeq(dds)
resultsNames(dds)

res <- results(dds, name="tissue_type_Tumor_vs_Normal")
library(EnhancedVolcano)

volcano_obj <- data.frame(
  "gene" = se@rowRanges@elementMetadata@listData[["gene_name"]],
  "logFC" = res@listData[["log2FoldChange"]],
  "P.Val" = res@listData[["padj"]]
)

EnhancedVolcano(volcano_obj,
                lab = volcano_obj$gene,
                x = 'logFC',
                pCutoff = 0.05,
                y = 'P.Val',
                drawConnectors = TRUE,
                legendPosition = 'right',
                labSize = 3.0,
                FCcutoff = 0.25,
                legendLabSize = 12,
                title = "Volcano Plot of LUAD",
                subtitle = "p-value cutoff = 0.05, log fold change cutoff=1",
                ylim = c(0,3)
)




