# ============================================================
# Multi-protein correlation analysis with BH multiple testing correction
# ============================================================

# 1. Load packages
library(tidyverse)
library(ggplot2)
library(corrplot)
library(psych)
library(ggpubr)
library(gridExtra)

# ============================================================
# 2. Define protein list and corresponding data files
# ============================================================

proteins <- c("HPD", "GLOD4", "GLO1", "MCEE")

# Store all results
all_results <- data.frame()

# ============================================================
# 3. Loop through each protein
# ============================================================

for (protein in proteins) {
  
  # Read data
  data_file <- paste0(protein, ".csv")
  data_raw <- read.csv(data_file, header = TRUE, stringsAsFactors = FALSE)
  
  # Data cleaning
  data_clean <- data_raw %>%
    select(Species, Length, Adult_survival, Max_age, Clutch_size, Egg_mass, Body_mass, ESI) %>%
    mutate(
      Length = as.numeric(as.character(Length)),
      Adult_survival = as.numeric(as.character(Adult_survival)),
      Max_age = as.numeric(as.character(Max_age)),
      Clutch_size = as.numeric(as.character(Clutch_size)),
      Egg_mass = as.numeric(as.character(Egg_mass)),
      Body_mass = as.numeric(as.character(Body_mass)),
      ESI = as.numeric(as.character(ESI))
    )
  
  # Select numeric variables
  data_numeric <- data_clean %>%
    select(Length, Adult_survival, Max_age, Clutch_size, Egg_mass, Body_mass, ESI) %>%
    mutate_all(as.numeric)
  
  # Remove variables with >50% missing values (retain Length and at least some traits)
  missing_prop <- colMeans(is.na(data_numeric))
  valid_vars <- names(missing_prop[missing_prop < 0.5])
  
  # Ensure Length is always in valid_vars
  if (!("Length" %in% valid_vars)) {
    valid_vars <- c("Length", valid_vars)
  }
  
  valid_data <- data_numeric[, valid_vars, drop = FALSE]
  
  # Get trait list (excluding Length)
  traits <- setdiff(valid_vars, "Length")
  
  # Perform correlation analysis for each trait
  for (trait in traits) {
    x <- valid_data$Length
    y <- valid_data[[trait]]
    
    # Remove missing values
    complete_idx <- complete.cases(x, y)
    x_clean <- x[complete_idx]
    y_clean <- y[complete_idx]
    n <- length(x_clean)
    
    if (n > 3) {
      test <- cor.test(x_clean, y_clean, method = "pearson")
      
      # Store results
      all_results <- rbind(all_results, data.frame(
        Protein = protein,
        Trait = trait,
        n = n,
        r = test$estimate,
        P_raw = test$p.value,
        CI_lower = test$conf.int[1],
        CI_upper = test$conf.int[2],
        stringsAsFactors = FALSE
      ))
    }
  }
}

# ============================================================
# 4. Multiple testing correction (Benjamini-Hochberg)
# ============================================================

all_results$P_adjusted <- p.adjust(all_results$P_raw, method = "BH")
all_results$Significant <- ifelse(all_results$P_adjusted < 0.05, "Yes", "No")

# Sort by protein and trait
all_results <- all_results %>%
  arrange(Protein, Trait)

# ============================================================
# 5. View significant results
# ============================================================

significant_results <- all_results[all_results$Significant == "Yes", ]
print(significant_results)

# ============================================================
# 6. Save results
# ============================================================

# Save full results (all 24 tests)
write.csv(all_results, "correlation_results_all_with_BH.csv", row.names = FALSE)

# Save significant results
if (nrow(significant_results) > 0) {
  write.csv(significant_results, "significant_results_after_BH.csv", row.names = FALSE)
}

# ============================================================
# 7. Print summary
# ============================================================

cat("\n========================================\n")
cat("Statistical summary:\n")
cat("Total tests:", nrow(all_results), "\n")
cat("Significant after BH correction:", nrow(significant_results), "\n")
cat("Significant results:\n")
print(significant_results)

# ============================================================
# 8. Generate publication-ready table (with significance markers)
# ============================================================

# Wide format table (protein × trait)
results_wide <- all_results %>%
  select(Protein, Trait, r, P_adjusted, n) %>%
  mutate(
    r_formatted = sprintf("%.2f", r),
    P_formatted = ifelse(P_adjusted < 0.001, "***",
                         ifelse(P_adjusted < 0.01, "**",
                                ifelse(P_adjusted < 0.05, "*", "ns"))),
    display = paste0(r_formatted, P_formatted)
  ) %>%
  select(Protein, Trait, display, n)

# Print table
print(results_wide)

# Save as CSV
write.csv(results_wide, "correlation_summary_table.csv", row.names = FALSE)