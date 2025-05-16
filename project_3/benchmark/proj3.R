library(tidyverse)
library(ggsignif)
library(broom)
library(ggeffects)

##### Functions #####
wilcox_boxplot_function <- function(vec1,
                                    vec2,
                                    group1_name = "Group1",
                                    group2_name = "Group2",
                                    alternative_hypothesis = "greater",
                                    y_axis_label = "Value",
                                    plot_title = "Boxplot with Wilcoxon Test") {
  
  # --- Perform one-tailed Wilcoxon rank-sum test ---
  alternative_hypothesis = ifelse(mean(vec1)>mean(vec2),"greater","less")
  
  fold_change = log2(mean(vec1)/mean(vec2))
  wilcox_result <- stats::wilcox.test(vec1, vec2,
                                      alternative = alternative_hypothesis,
                                      exact = NULL) # Let R decide on exact vs approximate
  
  p_value <- wilcox_result$p.value
  
  # --- Determine significance annotation ---
  significance_label <- dplyr::case_when(
    p_value <= 0.001 ~ "***",
    p_value <= 0.01  ~ "**",
    p_value <= 0.05  ~ "*",
    TRUE             ~ "ns" # Not significant
  )
  
  message(paste("P-value:", format(p_value, digits=4), "Significance:", significance_label))
  
  
  # --- Prepare data for ggplot ---
  data_long <- data.frame(
    Value = c(vec1, vec2),
    Group = factor(c(rep(group1_name, length(vec1)),
                     rep(group2_name, length(vec2))),
                   levels = c(group1_name, group2_name)) # Ensure correct order
  )
  
  # --- Create the boxplot using ggplot2 ---
  # Calculate y position for the significance bracket
  # This aims to place the bracket above the taller boxplot + some padding
  max_val_group1 <- if(length(vec1) > 0) max(fivenum(vec1)[5], quantile(vec1, 0.75) + 1.5 * IQR(vec1), na.rm = TRUE) else NA
  max_val_group2 <- if(length(vec2) > 0) max(fivenum(vec2)[5], quantile(vec2, 0.75) + 1.5 * IQR(vec2), na.rm = TRUE) else NA
  
  # Consider all data points for y_position calculation if outliers are present
  overall_max_val <- max(c(vec1, vec2), na.rm = TRUE)
  
  # Use a slightly more robust way to determine the top of the whiskers or outliers
  y_max_for_plot <- max(
    boxplot.stats(vec1)$stats[5], 
    boxplot.stats(vec2)$stats[5],
    if(any(boxplot.stats(vec1)$out)) max(boxplot.stats(vec1)$out) else -Inf,
    if(any(boxplot.stats(vec2)$out)) max(boxplot.stats(vec2)$out) else -Inf,
    na.rm = TRUE
  )
  
  # If all values are NA or infinite, y_max_for_plot could be -Inf.
  # Fallback to overall_max_val if calculation results in non-finite.
  if (!is.finite(y_max_for_plot)) {
    y_max_for_plot <- overall_max_val
  }
  
  # Add padding for the significance annotation
  y_position_signif <- y_max_for_plot + (0.1 * abs(y_max_for_plot))
  
  # If y_max_for_plot is 0 or very small, ensure y_position_signif is reasonably above
  if (y_max_for_plot == 0) {
    y_position_signif <- 0.1 * max(abs(data_long$Value), na.rm=TRUE)
    if (!is.finite(y_position_signif) || y_position_signif == 0) y_position_signif <- 1 # Default if all else fails
  } else if (abs(y_max_for_plot) < 1) { # For small values, a relative increase might be too small
    y_position_signif <- y_max_for_plot + 0.1 * max(abs(data_long$Value), na.rm=TRUE)
    if (!is.finite(y_position_signif)) y_position_signif <- y_max_for_plot + 0.2
  }
  
  
  plot_gg <- ggplot(data_long, aes(x = Group, y = Value, fill = Group)) +
    geom_boxplot(alpha = 0.7, outlier.shape = 21, outlier.colour = "grey30") +
    geom_jitter(width = 0.15, alpha = 0.2, shape = 16, size=1.5) + # Add individual points
    labs(title = paste(plot_title,"\nlog2FC:",format(fold_change,digits = 4)),
         x = "Group",
         y = y_axis_label) +
    theme_classic(base_size = 10) +
    theme(legend.position = "none",
          plot.title = element_text(hjust = 0.5, face = "bold"),
          axis.text = element_text(color = "black"),
          axis.title = element_text(face = "bold")) +
    scale_fill_brewer(palette = "Pastel1") # Example color palette
  
  # Add significance annotation only if not "ns" or if you always want the bar
  if (significance_label != "ns" || TRUE) { # Set to TRUE to always draw bar, even for "ns"
    plot_gg <- plot_gg +
      geom_signif(comparisons = list(c(group1_name, group2_name)),
                  annotations = significance_label,
                  y_position = y_position_signif,
                  map_signif_level = FALSE, # We are providing the label directly
                  tip_length = 0.02,
                  textsize = 5,
                  test = NULL, # We are not asking ggsignif to run a test
                  vjust = -0.2) # Adjust vertical position of the text
  }
  
  # Adjust y-axis limits to ensure the annotation is visible
  # Calculate a buffer for the y-axis limit based on the annotation position
  y_axis_upper_limit <- y_position_signif + (0.1 * abs(y_position_signif))
  # if y_position_signif is 0, use a small absolute buffer
  if (y_position_signif == 0) y_axis_upper_limit <- 0.2 * max(abs(data_long$Value), na.rm=TRUE)
  if (!is.finite(y_axis_upper_limit) || y_axis_upper_limit <= y_position_signif) {
    y_axis_upper_limit <- y_position_signif + 0.5 # Fallback buffer
  }
  
  current_y_limits <- layer_scales(plot_gg)$y$get_limits()
  plot_gg <- plot_gg + coord_cartesian(ylim = c(current_y_limits[1], max(current_y_limits[2], y_axis_upper_limit)))
  
  
  return(list(box_plot = plot_gg,
              hypo = alternative_hypothesis,
              p_value = p_value))
}

##### Analysis #####
###### Load Data ######
setwd("/Users/haowen/Documents/UCSD_course/BENG285")
# load mutation signatures
data_path = "/Users/haowen/Documents/GitHub/BENG285_team3/project_3/SigProf_Suggested_Solution/SBS96_De-Novo_Solution"
# list all folders
folders = list.dirs(data_path, full.names = TRUE, recursive = FALSE)
# get the path of signatures files like "*_Signatures.txt"
signature_files = list.files(folders, pattern = "*_Signatures.txt", full.names = TRUE, recursive = T)
# load all signatures
signature_list = lapply(signature_files, function(x) {
  read.table(x, header = TRUE, sep = "\t", row.names = 1)
})
merge_signature <- do.call(cbind,signature_list)

bootstrap_signature <- read.table("W_optim_k_4.tsv", header = TRUE, sep = "\t", row.names = 1)
bootstrap_signature <- bootstrap_signature/matrix(colSums(bootstrap_signature),
                                                  nrow = nrow(bootstrap_signature), 
                                                  ncol = ncol(bootstrap_signature),
                                                  byrow = T)
bootstrap_signature <- as.matrix(bootstrap_signature)
merge_signature <- as.matrix(merge_signature)
# correlation matrix
hm_mat <- crossprod(bootstrap_signature, merge_signature) / 
  (sqrt(colSums(bootstrap_signature^2)) %o% sqrt(colSums(merge_signature^2)))
library(ComplexHeatmap)
library(circlize)
Heatmap(hm_mat,cluster_columns = F, cluster_rows = F,
        col = colorRamp2(seq(0, 1,length.out=5), viridis::viridis(5)),
        rect_gp = gpar(col = "white", lwd = 1),column_title = "test")

## Two activities results
# get the path of signatures files like "*_Activities.txt"
activity_files = list.files(folders, pattern = "*_Activities_refit.txt", full.names = TRUE, recursive = T)
# load all signatures
activity_list = lapply(activity_files, function(x) {
  read.table(x, header = TRUE, sep = "\t", row.names = 1) #%>% 
  #`colnames<-`(paste0(colnames(.),"_", str_remove_all(str_extract(x,"_S[0-9]*_"),pattern = "_")))
})
merge_activity <- do.call(cbind,activity_list)

# from Aleysha
merge_activity <- read.table("H_optim_k_4.tsv", header = TRUE, sep = "\t", row.names = 1) 
# check the row names of activity_files are same
for (i in 1:length(activity_list)) {
  if (all(rownames(activity_list[[i]]) == rownames(activity_list[[1]])) == FALSE) {
    print(paste("The row names of activity file", i, "are not the same as the first one"))
  }
}
# extract the row names (samples) of the first activity file
rownames(merge_activity) <- str_sub(rownames(merge_activity),end = 12)
samples = rownames(merge_activity)
# check the row names of signature_files are same# load metadata
metadata = read.table("~/Documents/UCSD_course/BENG285/TCGA.LUAD.metadata.txt", header = TRUE, sep = "\t")
# subset the metadata to only include the samples in the activity files
metadata = metadata[metadata$patient_id %in% samples, ]
# load clinical data
#clinical = read_delim("~/Documents/UCSD_course/BENG285/clinical.tsv")
#clinical = clinical[clinical$cases.submitter_id %in% samples, ]

#pathology = read_delim("~/Documents/UCSD_course/BENG285/pathology_detail.tsv")
#pathology = pathology[pathology$cases.submitter_id %in% samples, ]

exposure = read_delim("~/Documents/UCSD_course/BENG285/exposure.tsv")
exposure = exposure[exposure$cases.submitter_id %in% samples, ]
# left join the metadata and clinical data
merge_metadata = left_join(metadata, exposure, 
                     by = c("patient_id" = "cases.submitter_id"), 
                     relationship = "one-to-many", multiple = "first")
# remove columns with only 1 unique value
merge_metadata = merge_metadata[, sapply(merge_metadata, function(x) length(unique(x)) > 1)]

merge_metadata$smoke_status <- sapply(merge_metadata$exposures.tobacco_smoking_status,
                                      \(s){
                                        if(s %in% c("Not Reported", "Unknown")){
                                          "NA"
                                        }else if(s == "Lifelong Non-Smoker"){
                                          "Non-smoker"
                                        }else{
                                          "Ever-smoker"
                                        }
                                      })
# same order
merge_activity <- merge_activity[merge_metadata$patient_id,]

nosmoke_idx <- merge_metadata$smoke_status == "Non-smoker"
smoke_idx <- merge_metadata$smoke_status == "Ever-smoker"
###### Box Plot ######
# based on activity, find signature significant different in different group
diff_sig_df <- data.frame(sig_id = colnames(merge_activity),
                          p_val = 1,
                          hypo = "NA",
                          stringsAsFactors = FALSE)
boxplot_ls <- list()
for(i in 1:nrow(diff_sig_df)){
  ctrl_vec <- log10(1+merge_activity[nosmoke_idx,diff_sig_df$sig_id[i]])
  smoke_vec <- log10(1+merge_activity[smoke_idx,diff_sig_df$sig_id[i]])
  
  tmp <- wilcox_boxplot_function(smoke_vec,ctrl_vec,
                                 group1_name = "Ever-smoker",
                                 group2_name = "Non-smoker",
                                 plot_title = diff_sig_df$sig_id[i])
  diff_sig_df[i,"p_value"] = tmp$p_value
  diff_sig_df[i,"hypo"] = tmp$hypo
  
  if(tmp$p_value < 0.05){
    boxplot_ls[[diff_sig_df$sig_id[i]]] <- tmp$box_plot
  }
  
}

ggpubr::ggarrange(plotlist = boxplot_ls,ncol = 5, nrow = 1)

###### Corr w/ Age ######
scatter_plot_ls <- list()
for(i in 1:nrow(diff_sig_df)){
  tmp_df <- data.frame(ma = log10(1+merge_activity[,diff_sig_df$sig_id[i]]),
                       age = merge_metadata$age_at_initial_pathologic_diagnosis,
                       smoke = merge_metadata$smoke_status) %>% subset(smoke!="NA"&!is.na(age))
  tmp_df$smoke <- as.factor(tmp_df$smoke)
  glm_model <- glm(ma ~ age + smoke,
                   data = tmp_df,
                   family = gaussian()) 
  
  age_coeffs <- tidy(glm_model) %>% filter(term == "age")
  age_coefficient_value <- age_coeffs$estimate
  age_p_value <- age_coeffs$p.value
  
  age_effect_preds <- ggpredict(glm_model, terms = "age")
  
  plot_title_text <- paste0(diff_sig_df$sig_id[i],
                            " - Age Coeff: ", format(age_coefficient_value, digits = 3),
                            "; P-value: ", format(age_p_value, digits = 3))
  
  scatter_plot_ls[[diff_sig_df$sig_id[i]]] <- ggplot() +
    # Scatter plot of raw data, color by smoking status for context
    geom_point(data = tmp_df, aes(x = age, y = ma, color = smoke), alpha = 0.7) +
    # Add the fitting curve for the age effect (adjusted)
    geom_line(data = age_effect_preds, aes(x = x, y = predicted), color = "blue", linewidth = 1.2) +
    geom_ribbon(data = age_effect_preds, aes(x = x, ymin = conf.low, ymax = conf.high),
                alpha = 0.2, fill = "blue") +
    scale_color_manual(values = c("Ever-smoker" = "red", "Non-smoker" = "cornflowerblue"),
                       name = "Smoking Status") +
    labs(title = plot_title_text,
         x = "Age (Years)",
         y = "Log10 Mutation Activity") +
    theme_minimal(base_size = 10) +
    theme(plot.title = element_text(hjust = 0.5, size = 12),
          legend.position = "top")
  
}

ggpubr::ggarrange(plotlist = scatter_plot_ls[-2],ncol = 2, nrow = 2)



