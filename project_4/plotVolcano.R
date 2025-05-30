library(EnhancedVolcano)
data <- read_csv("./deseq2_results_full.csv")
volcano_obj <- data.frame(
  "HGNC_Symbol" = data$...1,
  "logFC" = data$log2FoldChange,
  "adj.P.Val" = data$padj
)

# plot with package EnhancedVolcano
EnhancedVolcano(volcano_obj,
                lab = volcano_obj$HGNC_Symbol,
                x = 'logFC',
                pCutoff = 0.05,
                y = 'adj.P.Val',
                drawConnectors = TRUE,
                legendPosition = 'left',
                labSize = 3.0,
                legendLabSize = 12,
                title = "Volcano Plot",
                subtitle = "adj.p-value cutoff = 0.05, log fold change cutoff=1",
                ylim = c(0,3)
)
