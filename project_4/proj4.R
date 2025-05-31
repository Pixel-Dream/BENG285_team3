library(survival)
library(survminer) 

# 2. Simulate Example Data
# !!! IMPORTANT: Replace this section with your actual data loading and preparation !!!
set.seed(123) # for reproducibility

df <- read_csv("Documents/UCSD_course/BENG285/TCGA_LUAD_final_for_OS_model.csv") %>%
  column_to_rownames("patient_id")

df_mut <- read_csv("Documents/UCSD_course/BENG285/KM_mutation_top10.csv") %>%
  column_to_rownames("patient_id") %>%
  `colnames<-`(paste0(colnames(.),"_M"))

df_de <- read_csv("Documents/UCSD_course/BENG285/KM_DEG_top10_expression.csv") %>%
  group_by(patient_id) %>%
  summarise(across(all_of(colnames(.)[-1]), mean, .names = "{.col}")) %>%
  ungroup() %>%
  column_to_rownames("patient_id") %>%
  `colnames<-`(paste0(colnames(.),"_E"))

df <- tibble(df,df_mut[row.names(df),],df_de[row.names(df),-c(1,2)])


# Ensure your categorical variables are factors. 'stage' might be ordered if appropriate.
# df$gender <- factor(df$gender)
# df$race <- factor(df$race)
# df$stage <- factor(df$stage, ordered = TRUE, levels = c("I", "II", "III", "IV")) # If stage is ordinal

# 3. Define variable groups
prognostic_covariates <- c("age", "gender", "race", "stage")

valid_colnames <- make.names(colnames(df))
colnames(df) <- valid_colnames


#all_gene_predictors <- colnames(df)[9:28]
all_gene_predictors <- colnames(df)[c(9,11:16,18:21,23:28)]
all_model_predictors <- c(prognostic_covariates, all_gene_predictors)

# Check for missing data and handle it (e.g., imputation, removal)
# For this example, we'll proceed assuming data is complete or handled.
# print(sapply(df[, c("OS.time", "OS", all_model_predictors)], function(x) sum(is.na(x))))
# df_complete <- na.omit(df[, c("OS.time", "OS", all_model_predictors)]) # Listwise deletion example

# Create the survival object and formula string
formula_string <- paste("Surv(OS.time, OS) ~ ", paste(all_model_predictors, collapse = " + "))
cox_formula <- as.formula(formula_string)

cat("## Cox Model Formula:\n")
print(cox_formula)

# 5. Fit the Cox Proportional Hazards Model
# Using 'efron' for tie handling is often a good default.
# It's crucial to ensure factors have more than one level present in the data.
# For example, if after filtering, a race category has 0 patients, it can cause errors.
cox_model <- tryCatch({
  coxph(cox_formula, data = as.data.frame(df), ties = "efron",x = T)
}, error = function(e) {
  message("Error fitting Cox model: ", e$message)
  # You might want to investigate problematic variables if an error occurs:
  # for(p in all_model_predictors) {
  #   if(is.factor(df[[p]])) { message(p, ": ", paste(levels(df[[p]]), collapse=", "), " (", length(unique(df[[p]])), " unique levels in data)"); print(table(df[[p]])) }
  #   else { message(p, ": ", length(unique(df[[p]])), " unique values") }
  # }
  return(NULL)
})

cat("\n## Cox Model Summary:\n")
model_summary <- summary(cox_model)
print(model_summary)

# 7. Address Multiple Hypothesis Testing for Gene Predictors
cat("\n## Multiple Hypothesis Testing for Gene Predictors:\n")
coef_summary <- model_summary$coefficients

# Identify gene predictors present in the model output
# Factor variables will have names like 'raceBlack', 'raceAsian' etc.
# Numeric variables (like gene expressions or 0/1 mutations) will have their original names.

# We are interested in the p-values for the 20 gene-related variables
# Ensure the names match exactly how they appear in `rownames(coef_summary)`

gene_predictors_in_model <- intersect(all_gene_predictors, rownames(coef_summary))

if (length(gene_predictors_in_model) > 0) {
  gene_p_values <- coef_summary[gene_predictors_in_model, "Pr(>|z|)"] # Last column of coefficients table
  
  # Bonferroni correction (more conservative)
  p_values_bonferroni <- p.adjust(gene_p_values, method = "bonferroni")
  
  # Benjamini-Hochberg (FDR) correction (often preferred)
  p_values_fdr <- p.adjust(gene_p_values, method = "BH")
  
  adjusted_p_results <- data.frame(
    Predictor = gene_predictors_in_model,
    Hazard_Ratio = exp(coef_summary[gene_predictors_in_model, "coef"]),
    Original_P_Value = gene_p_values,
    Bonferroni_Adjusted_P = p_values_bonferroni,
    FDR_Adjusted_P = p_values_fdr,
    stringsAsFactors = FALSE
  )
  rownames(adjusted_p_results) <- NULL # Clean row names
  
  cat("### Adjusted P-values for Gene Predictors:\n")
  print(adjusted_p_results)
  
  alpha <- 0.05 # Significance level
  significant_genes_fdr <- adjusted_p_results[adjusted_p_results$FDR_Adjusted_P < alpha, ]
  cat(paste0("\n### Gene Predictors Significant at FDR < ", alpha, ":\n"))
  if(nrow(significant_genes_fdr) > 0){
    print(significant_genes_fdr)
  } else {
    cat("No gene predictors were significant after FDR correction.\n")
  }
  
} else {
  cat("No gene-related predictors found directly in the model coefficient matrix for p-value adjustment.\n")
  cat("This could be due to variable removal, naming mismatches, or other model fitting issues.\n")
  cat("Model coefficient row names:\n")
  print(rownames(coef_summary))
}



model_summary <- summary(cox_model)
coef_table <- model_summary$coefficients
conf_int_table <- model_summary$conf.int

# Define groups and their members (match names from coef_table rownames)
# Rownames for factors will be like 'genderMale', 'raceBlack', 'stageII'
# Numeric variables like 'age', 'MROH9' will have their direct names.

# --- Define groups and map variable names ---
# User's desired grouping structure. We need to find these in `rownames(coef_table)`.
# Note: Reference categories (e.g., genderFemale, raceWhite, stageI) won't appear in coef_table.

# Helper to find model coefficient names (case-insensitive for flexibility with user input)
find_coef_name <- function(pattern, all_names) {
  match <- grep(paste0("^", pattern, "$"), all_names, ignore.case = TRUE, value = TRUE)
  if (length(match) == 1) return(match)
  # For factors, pattern might be 'genderMALE' and name is 'genderMale'
  match_factor <- grep(pattern, all_names, ignore.case = TRUE, value = TRUE)
  if (length(match_factor) == 1) return(match_factor)
  else if(length(match_factor) > 1) return(match_factor[1])
  return(NA_character_) # Return NA if not found or ambiguous
}

all_model_coef_names <- rownames(coef_table)

# Group definitions
# The names used here (e.g. "genderMale") MUST match the output from the cox model summary
# For factors, this is typically variable_name + level_name (without spaces)
group_defs <- list(
  "Age" = find_coef_name("age", all_model_coef_names),
  "Gender" = find_coef_name("genderMale", all_model_coef_names), # Assuming Female is reference
  "Race" = c(find_coef_name("raceASIAN", all_model_coef_names), 
             find_coef_name("raceBLACK OR AFRICAN AMERICAN", all_model_coef_names),
             find_coef_name("raceWHITE", all_model_coef_names)
  ),
  "Stage" = c(find_coef_name("Stage II", all_model_coef_names),
              find_coef_name("Stage III", all_model_coef_names),
              find_coef_name("Stage IV", all_model_coef_names)
              # Note: stageI is reference
  ),
  "Mutation" = sapply(c("MROH9","BSN","MAP3K19","ABCA12","ZNF521","SMARCA4","STK11","PLXNA2"), 
                      find_coef_name, all_names = all_model_coef_names, USE.NAMES = FALSE),
  "DEG Expression" = sapply(c("RSPO1","DIO1","TRIM50","IGLON5","GAL","PKHD1L1","HSPA6","CREG2","PVALB"),
                            find_coef_name, all_names = all_model_coef_names, USE.NAMES = FALSE)
)

# Clean up NAs from group_defs (variables not found in model output)
group_defs <- lapply(group_defs, function(x) x[!is.na(x)])
group_defs <- group_defs[sapply(group_defs, length) > 0 | names(group_defs) %in% c("Age", "Gender")] # Keep headers even if var not found for age/gender for structure

# --- Prepare data for forestplot ---
plot_elements_list <- list()
is_summary_list <- list()
text_color_list <- list()

# Significance level
alpha <- 0.05

for (group_name in names(group_defs)) {
  # Add group header
  plot_elements_list[[length(plot_elements_list) + 1]] <- list(
    name = group_name, hr_text = NA, p_val_text = NA, 
    mean = NA, lower = NA, upper = NA, p_value = NA
  )
  is_summary_list[[length(is_summary_list) + 1]] <- TRUE
  text_color_list[[length(text_color_list) + 1]] <- "black" # Header color
  
  # Add variables in the group
  vars_in_group <- group_defs[[group_name]]
  if (length(vars_in_group) > 0 && !all(is.na(vars_in_group))) {
    for (var_name in vars_in_group) {
      if (var_name %in% rownames(coef_table)) {
        p_value <- coef_table[var_name, "Pr(>|z|)"]
        plot_elements_list[[length(plot_elements_list) + 1]] <- list(
          name = paste("   ", var_name), # Indent variable names
          hr_text = sprintf("%.2f (%.2f-%.2f)", 
                            conf_int_table[var_name, "exp(coef)"], 
                            conf_int_table[var_name, "lower .95"], 
                            conf_int_table[var_name, "upper .95"]),
          p_val_text = format.pval(p_value, digits = 2, eps = 0.001),
          mean = conf_int_table[var_name, "exp(coef)"],
          lower = conf_int_table[var_name, "lower .95"],
          upper = conf_int_table[var_name, "upper .95"],
          p_value = p_value
        )
        is_summary_list[[length(is_summary_list) + 1]] <- FALSE
        text_color_list[[length(text_color_list) + 1]] <- ifelse(p_value < alpha, "red", "black")
      }
    }
  }
}

# Convert list of lists to a data frame for easier extraction
plot_df <- do.call(rbind.data.frame, lapply(plot_elements_list, function(x) as.data.frame(t(unlist(x)), stringsAsFactors=FALSE)))
# Ensure numeric columns are numeric
numeric_cols <- c("mean", "lower", "upper", "p_value")
for(col in numeric_cols) {
  plot_df[[col]] <- as.numeric(plot_df[[col]])
}


# Prepare tabletext for forestplot
tabletext_matrix <- cbind(
  c("Predictor", plot_df$name),
  c("Hazard Ratio (95% CI)", plot_df$hr_text),
  c("P-value", plot_df$p_val_text)
)
tabletext_matrix[is.na(tabletext_matrix)] <- "" # Replace NA with empty string for display

# Prepare plot data (mean, lower, upper) - first row is for header, so NA
plot_mean <- c(NA, plot_df$mean)
plot_lower <- c(NA, plot_df$lower)
plot_upper <- c(NA, plot_df$upper)

# is.summary vector - first row is header, then based on list
is_summary_vector <- c(TRUE, unlist(is_summary_list)) # TRUE for main header, then from list

# Text colors vector - first row is header, then from list
text_colors_vector <- c("black", unlist(text_color_list)) # Main header black

# --- Data Sanity Checks and Clipping (same as before) ---
min_hr_plot <- 0.05 # Adjusted for potentially more variability
max_hr_plot <- 10  # Adjusted

# Apply to plot_mean, plot_lower, plot_upper, skipping NAs (headers)
non_na_indices <- !is.na(plot_mean)
plot_mean_valid <- plot_mean[non_na_indices]
plot_lower_valid <- plot_lower[non_na_indices]
plot_upper_valid <- plot_upper[non_na_indices]

# Replace NA/NaN/Inf for plotting values
plot_mean_valid[is.na(plot_mean_valid) | is.nan(plot_mean_valid) | is.infinite(plot_mean_valid)] <- 1
plot_lower_valid[is.na(plot_lower_valid) | is.nan(plot_lower_valid) | plot_lower_valid <= 0 | is.infinite(plot_lower_valid)] <- min_hr_plot
plot_upper_valid[is.na(plot_upper_valid) | is.nan(plot_upper_valid) | plot_upper_valid <= 0 | is.infinite(plot_upper_valid)] <- max_hr_plot

# Ensure mean is within bounds
plot_mean_valid <- pmax(min_hr_plot, pmin(max_hr_plot, plot_mean_valid))

# Ensure lower < mean < upper after clipping
plot_lower_valid <- pmin(plot_lower_valid, plot_mean_valid)
plot_upper_valid <- pmax(plot_upper_valid, plot_mean_valid)

# Put cleaned values back, handling NAs for headers
plot_mean[non_na_indices] <- plot_mean_valid
plot_lower[non_na_indices] <- plot_lower_valid
plot_upper[non_na_indices] <- plot_upper_valid

cat("\n### Final plot data for forestplot (mean, lower, upper - first row is NA for header):\n")
print(data.frame(mean=plot_mean, lower=plot_lower, upper=plot_upper))
cat("\n### Table text:\n")
print(tabletext_matrix)
cat("\n### Is Summary vector:\n")
print(is_summary_vector)
cat("\n### Text Colors vector:\n")
print(text_colors_vector)


# --- Perform forest plot ---
if (nrow(plot_df) > 0) {
  plot_clip_range <- c(min_hr_plot, max_hr_plot)
  
  # Define ticks, ensuring they are within the clip range and sensible
  # Example ticks:
  # potential_ticks <- c(0.1, 0.25, 0.5, 0.75, 1, 1.5, 2, 5, 10)
  # xticks_filtered <- potential_ticks[potential_ticks >= plot_clip_range[1] & potential_ticks <= plot_clip_range[2]]
  # if(length(xticks_filtered) < 2 || any(is.na(xticks_filtered))) xticks_filtered <- NULL 
  # A simpler auto approach for ticks based on clipped range:
  xticks_filtered <- pretty(plot_clip_range, n=5) # Generate ~5 nice ticks within the range
  xticks_filtered <- xticks_filtered[xticks_filtered > 0] # Ensure ticks are positive for log scale
  if(length(xticks_filtered) < 2) xticks_filtered <- NULL
  
  require(forestplot)
  forestplot::forestplot(
    labeltext = tabletext_matrix,
    mean = plot_mean,
    lower = plot_lower,
    upper = plot_upper,
    is.summary = is_summary_vector,
    title = "Forest Plot of Hazard Ratios by Group",
    txt_gp = fpTxtGp(label = grid::gpar(cex = 0.75), # Text for labels
                     ticks = grid::gpar(cex = 0.7),  # Tick mark numbers
                     xlab  = grid::gpar(cex = 0.9),  # X-axis label
                     title = grid::gpar(cex = 1.1)), # Plot title
    col = fpColors(box = ifelse(text_colors_vector[-1] == "red" & !is_summary_vector[-1], "darkred", "royalblue"), # Highlight significant boxes
                   line = ifelse(text_colors_vector[-1] == "red" & !is_summary_vector[-1], "red", "darkblue"),   # Highlight significant lines
                   summary = "black", # Color for summary lines (group headers)
                   text = text_colors_vector), # Vector of colors for text
    boxsize = 0.2,
    graph.pos = ncol(tabletext_matrix), # Graph after the last text column
    hrzl_lines = gpar(col="#999999", lty=2), # Style for horizontal lines between vars
    # Add lines to separate major groups visually
    # This requires knowing the row indices of the group headers in the final plot structure
    # For now, using the default hrzl_lines. More complex lines would need manual index calculation.
    clip = plot_clip_range,
    xticks = xticks_filtered,
    mar = grid::unit(c(5,2,5,2), "mm"), # c(bottom, left, top, right) margins
    xlab = "Hazard Ratio (95% CI)"
  )
} else {
  cat("No data available to plot after grouping and filtering.\n")
}