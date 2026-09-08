# ============================================================
# Generate training dataset from HPD multiple sequence alignment
# Input:
#   - avian_hpd_data.csv: contains ID, Length, Sequence, VOC2_C_position, Class
#   - avian_hpd_alignment.fasta: multiple sequence alignment file (FASTA format)
#   - Reference sequence ID: used to define the coordinate system
# Output:
#   - hpd_training_data.csv: features for each alignment position
#   - Two PDF plots: global extension potential + branch comparison
# ============================================================

# Load required packages
library(tidyverse)
library(Biostrings)

# ============================================================
# 1. Set parameters
# ============================================================

# File paths
data_csv <- "avian_hpd_data.csv"
alignment_fasta <- "avian_hpd_alignment.fasta"

# Reference sequence ID (must exist in the data table)
ref_id <- "A0A6J0H6Y4"  # Replace with the actual reference sequence ID

# Output file
output_csv <- "hpd_training_data.csv"

# ============================================================
# 2. Read sequence information table
# ============================================================

species_data <- read.csv(data_csv, stringsAsFactors = FALSE)

# Check required columns
required_cols <- c("ID", "Length", "Sequence", "VOC2_C_position", "Class")
missing_cols <- required_cols[!required_cols %in% colnames(species_data)]
if (length(missing_cols) > 0) {
  stop("Missing columns in data table: ", paste(missing_cols, collapse = ", "))
}

# Check Class column for valid values
valid_classes <- c("main", "alternative")
invalid_classes <- setdiff(unique(species_data$Class), valid_classes)
if (length(invalid_classes) > 0) {
  warning("Unexpected values in Class column: ", paste(invalid_classes, collapse = ", "))
}

cat("Data summary:\n")
cat("  Total sequences: ", nrow(species_data), "\n")
cat("  main branch: ", sum(species_data$Class == "main"), "\n")
cat("  alternative branch: ", sum(species_data$Class == "alternative"), "\n")

# Check if reference sequence exists
if (!(ref_id %in% species_data$ID)) {
  stop("Reference sequence ID not found in data table")
}

# Extract reference sequence information
ref_info <- species_data %>% filter(ID == ref_id)
ref_length <- ref_info$Length
ref_voc2_c <- ref_info$VOC2_C_position
ref_class <- ref_info$Class

cat("Reference sequence information:\n")
cat("  ID: ", ref_id, "\n")
cat("  Branch: ", ref_class, "\n")
cat("  Length: ", ref_length, "\n")
cat("  VOC2_C position: ", ref_voc2_c, "\n")

# ============================================================
# 3. Read multiple sequence alignment
# ============================================================

alignment <- readAAStringSet(alignment_fasta)
seq_names <- names(alignment)

# Ensure sequence names match IDs in data table
if (!all(seq_names %in% species_data$ID)) {
  warning("Some sequence IDs not found in data table, keeping only matching sequences")
  alignment <- alignment[seq_names %in% species_data$ID]
  seq_names <- names(alignment)
}

# Convert to character matrix (rows = sequences, columns = positions)
alignment_matrix <- as.matrix(alignment)
n_seq <- nrow(alignment_matrix)
n_pos <- ncol(alignment_matrix)

cat("Alignment matrix: ", n_seq, " sequences, ", n_pos, " positions\n")

# ============================================================
# 4. Locate reference sequence row index
# ============================================================

ref_idx <- which(seq_names == ref_id)
if (length(ref_idx) == 0) {
  stop("Reference sequence ID not found in alignment file")
}

# Retrieve Class labels for each sequence
seq_classes <- species_data$Class[match(seq_names, species_data$ID)]
seq_classes <- as.character(seq_classes)

# Get indices for main and alternative branches
main_idx <- which(seq_classes == "main")
alt_idx <- which(seq_classes == "alternative")

# ============================================================
# 5. Calculate statistical features for each position (including branch-specific)
# ============================================================

results <- data.frame(
  reference_position = 1:n_pos,
  distance_from_VOC2_C = NA_integer_,
  # Global statistics
  presence_frequency = NA_real_,
  extension_potential = NA_real_,
  gap_frequency = NA_real_,
  # main branch statistics
  main_presence_frequency = NA_real_,
  main_extension_potential = NA_real_,
  # alternative branch statistics
  alt_presence_frequency = NA_real_,
  alt_extension_potential = NA_real_,
  # Sequence information
  reference_aa = NA_character_,
  most_frequent_aa = NA_character_,
  aa_conservation_score = NA_real_,
  aa_diversity = NA_real_,
  stringsAsFactors = FALSE
)

# Progress bar
pb <- txtProgressBar(min = 0, max = n_pos, style = 3)

for (i in 1:n_pos) {
  col_chars <- alignment_matrix[, i]
  
  # ----- Global statistics -----
  non_gap <- col_chars[col_chars != "-"]
  n_non_gap <- length(non_gap)
  presence_freq <- n_non_gap / n_seq
  gap_freq <- 1 - presence_freq
  extension_potential <- gap_freq
  
  # ----- main branch statistics -----
  if (length(main_idx) > 0) {
    main_chars <- col_chars[main_idx]
    main_non_gap <- main_chars[main_chars != "-"]
    main_presence <- length(main_non_gap) / length(main_idx)
  } else {
    main_presence <- NA_real_
  }
  main_extension <- 1 - main_presence
  
  # ----- alternative branch statistics -----
  if (length(alt_idx) > 0) {
    alt_chars <- col_chars[alt_idx]
    alt_non_gap <- alt_chars[alt_chars != "-"]
    alt_presence <- length(alt_non_gap) / length(alt_idx)
  } else {
    alt_presence <- NA_real_
  }
  alt_extension <- 1 - alt_presence
  
  # ----- Reference sequence amino acid -----
  ref_aa <- ifelse(col_chars[ref_idx] != "-", col_chars[ref_idx], NA_character_)
  
  # ----- Amino acid frequency statistics (global) -----
  if (n_non_gap > 0) {
    aa_counts <- table(non_gap)
    most_freq <- names(aa_counts)[which.max(aa_counts)]
    max_freq <- max(aa_counts) / n_non_gap
    freq_vector <- aa_counts / n_non_gap
    shannon <- -sum(freq_vector * log(freq_vector + 1e-10))
  } else {
    most_freq <- NA_character_
    max_freq <- NA_real_
    shannon <- NA_real_
  }
  
  # Store results
  results$distance_from_VOC2_C[i] <- i - ref_voc2_c
  results$presence_frequency[i] <- presence_freq
  results$extension_potential[i] <- extension_potential
  results$gap_frequency[i] <- gap_freq
  results$main_presence_frequency[i] <- main_presence
  results$main_extension_potential[i] <- main_extension
  results$alt_presence_frequency[i] <- alt_presence
  results$alt_extension_potential[i] <- alt_extension
  results$reference_aa[i] <- ref_aa
  results$most_frequent_aa[i] <- most_freq
  results$aa_conservation_score[i] <- max_freq
  results$aa_diversity[i] <- shannon
  
  setTxtProgressBar(pb, i)
}
close(pb)

cat("Statistical calculation completed\n")

# ============================================================
# 6. Add amino acid physicochemical properties
# ============================================================

aa_properties <- data.frame(
  aa = c("A","R","N","D","C","Q","E","G","H","I","L","K","M","F","P","S","T","W","Y","V"),
  group = c("Hydrophobic","Positive","Polar","Negative","Polar","Polar","Negative","Polar",
            "Positive","Hydrophobic","Hydrophobic","Positive","Hydrophobic","Hydrophobic",
            "Hydrophobic","Polar","Polar","Hydrophobic","Polar","Hydrophobic"),
  hydrophobicity = c(1.8, -4.5, -3.5, -3.5, 2.5, -3.5, -3.5, -0.4, -3.2, 4.5, 3.8, -3.9,
                     1.9, 2.8, -1.6, -0.8, -0.7, -0.9, -1.3, 4.2),
  charge = c("neutral","positive","neutral","negative","neutral","neutral","negative","neutral",
             "positive","neutral","neutral","positive","neutral","neutral","neutral","neutral",
             "neutral","neutral","neutral","neutral"),
  molecular_weight = c(89.09, 174.20, 132.12, 133.10, 121.16, 146.15, 147.13, 75.07,
                       155.16, 131.17, 131.17, 146.19, 149.21, 165.19, 115.13, 105.09,
                       119.12, 204.23, 181.19, 117.15)
)

results <- merge(results, aa_properties, by.x = "reference_aa", by.y = "aa", all.x = TRUE)
colnames(results)[colnames(results) == "group"] <- "aa_group"
colnames(results)[colnames(results) == "hydrophobicity"] <- "aa_hydrophobicity"
colnames(results)[colnames(results) == "charge"] <- "aa_charge"
colnames(results)[colnames(results) == "molecular_weight"] <- "aa_molecular_weight"

# ============================================================
# 7. Add window-based features
# ============================================================

window_size <- 5
results$window_presence_mean <- NA_real_
results$window_gap_mean <- NA_real_

for (i in 1:n_pos) {
  start_idx <- max(1, i - window_size)
  end_idx <- min(n_pos, i + window_size)
  window_presence <- results$presence_frequency[start_idx:end_idx]
  results$window_presence_mean[i] <- mean(window_presence, na.rm = TRUE)
  results$window_gap_mean[i] <- mean(1 - window_presence, na.rm = TRUE)
}

# ============================================================
# 8. Generate final data and export
# ============================================================

results$site_id <- paste0("HPD_", results$reference_position)

final_columns <- c(
  "site_id",
  "reference_position",
  "distance_from_VOC2_C",
  # Global statistics
  "presence_frequency",
  "extension_potential",
  "gap_frequency",
  # main branch statistics
  "main_presence_frequency",
  "main_extension_potential",
  # alternative branch statistics
  "alt_presence_frequency",
  "alt_extension_potential",
  # Sequence features
  "reference_aa",
  "most_frequent_aa",
  "aa_conservation_score",
  "aa_diversity",
  "aa_group",
  "aa_hydrophobicity",
  "aa_charge",
  "aa_molecular_weight",
  # Window features
  "window_presence_mean",
  "window_gap_mean"
)

final_data <- results[, final_columns]

write.csv(final_data, output_csv, row.names = FALSE)

cat("\nTraining dataset generated:", output_csv, "\n")
cat("Total rows:", nrow(final_data), "\n")
cat("Reference sequence:", ref_id, "\n")
cat("Reference sequence branch:", ref_class, "\n")

# ============================================================
# 9. Plot: global extension potential
# ============================================================

library(ggplot2)

p1 <- ggplot(final_data, aes(x = distance_from_VOC2_C, y = extension_potential)) +
  geom_line(color = "steelblue", size = 1) +
  geom_point(alpha = 0.3, size = 0.5) +
  geom_hline(yintercept = 0.5, linetype = "dashed", color = "gray50", alpha = 0.5) +
  labs(title = "Global extension potential: all Aves species",
       x = "Distance from VOC2_C (negative = N-terminus, positive = C-terminus)",
       y = "Extension potential (1 - presence frequency)") +
  theme_minimal(base_size = 12)

ggsave("extension_potential_global.pdf", plot = p1, width = 10, height = 6)
cat("Saved: extension_potential_global.pdf\n")

# ============================================================
# 10. Plot: branch comparison of extension potential
# ============================================================

# Convert to long format for grouped plotting
branch_comparison <- final_data %>%
  select(distance_from_VOC2_C, 
         main_extension_potential, 
         alt_extension_potential) %>%
  pivot_longer(
    cols = c(main_extension_potential, alt_extension_potential),
    names_to = "branch",
    values_to = "extension_potential"
  ) %>%
  mutate(branch = ifelse(branch == "main_extension_potential", "Main", "Alternative"))

p2 <- ggplot(branch_comparison, aes(x = distance_from_VOC2_C, y = extension_potential, color = branch)) +
  geom_line(size = 1) +
  geom_hline(yintercept = 0.5, linetype = "dashed", color = "gray50", alpha = 0.5) +
  scale_color_manual(values = c("Main" = "#2E86C1", "Alternative" = "#E74C3C")) +
  labs(title = "Branch comparison: Main vs Alternative extension potential",
       x = "Distance from VOC2_C (negative = N-terminus, positive = C-terminus)",
       y = "Extension potential (1 - presence frequency)",
       color = "Branch") +
  theme_minimal(base_size = 12) +
  theme(legend.position = "top")

ggsave("extension_potential_branch_comparison.pdf", plot = p2, width = 10, height = 6)
cat("Saved: extension_potential_branch_comparison.pdf\n")

# ============================================================
# 11. Data summary
# ============================================================

cat("\nData summary:\n")
cat("  distance range: ", range(final_data$distance_from_VOC2_C, na.rm = TRUE), "\n")
cat("  Global extension potential range: ", range(final_data$extension_potential, na.rm = TRUE), "\n")
cat("  main branch extension potential range: ", range(final_data$main_extension_potential, na.rm = TRUE), "\n")
cat("  alternative branch extension potential range: ", range(final_data$alt_extension_potential, na.rm = TRUE), "\n")

cat("\n✅ Complete!\n")