# ============================================================
# HPD extension feature analysis visualization
# Input: hpd_training_data.csv, avian_hpd_data.csv, avian_hpd_alignment.fasta
# Output: 3 PDF figures (Fig1, Fig2, Fig3)
# ============================================================

# Load required packages
library(tidyverse)
library(Biostrings)
library(ggplot2)

# ============================================================
# 1. Load data
# ============================================================

# Position-level summary data
df_pos <- read.csv("hpd_training_data.csv", stringsAsFactors = FALSE)

# Species data (for length distribution)
species_info <- read.csv("avian_hpd_data.csv", stringsAsFactors = FALSE)

# Multiple sequence alignment (for computing N-terminal gap ratio)
alignment <- readAAStringSet("avian_hpd_alignment.fasta")
alignment_matrix <- as.matrix(alignment)
seq_names <- rownames(alignment_matrix)

# Ensure sequence order matches species_info
species_info <- species_info[match(seq_names, species_info$ID), ]

# Reference sequence ID and VOC2_C position
ref_id <- "A0A6J0H6Y4"  # Replace with your reference ID
ref_idx <- which(seq_names == ref_id)
ref_voc2_c <- species_info$VOC2_C_position[species_info$ID == ref_id]

# Distance vector (matches distance_from_VOC2_C in summary data)
distance_vec <- 1:ncol(alignment_matrix) - ref_voc2_c

# ============================================================
# Figure 1: Extension depth profile along the sequence (core figure)
# ============================================================

p1 <- ggplot(df_pos, aes(x = distance_from_VOC2_C, y = extension_potential)) +
  geom_line(color = "steelblue", size = 1.2) +
  geom_ribbon(aes(ymin = extension_potential - 0.02, ymax = extension_potential + 0.02),
              alpha = 0.2, fill = "steelblue") +
  geom_vline(xintercept = 0, linetype = "dashed", color = "red", alpha = 0.6) +
  annotate("text", x = -50, y = 0.9, label = "N-terminal extension", color = "blue", size = 5) +
  annotate("text", x = 50, y = 0.9, label = "C-terminal conserved core", color = "red", size = 5) +
  labs(title = "Extension depth profile along the sequence",
       x = "Distance from VOC2_C (negative = N-terminus, positive = C-terminus)",
       y = "Extension depth (gap frequency)") +
  theme_minimal(base_size = 14) +
  theme(plot.title = element_text(hjust = 0.5))

# ============================================================
# Figure 2: Feature comparison between conserved and extension regions (boxplot)
# ============================================================

# Define two regions: Conserved (>0) and Extension (<0)
df_pos <- df_pos %>%
  mutate(region = ifelse(distance_from_VOC2_C > 0, "Conserved", "Extension"))

# Select three features: conservation, diversity, hydrophobicity
features_long <- df_pos %>%
  select(region, aa_conservation_score, aa_diversity, aa_hydrophobicity) %>%
  pivot_longer(cols = -region, names_to = "feature", values_to = "value")

# Feature labels (English, for PDF display)
feature_labels <- c(
  "aa_conservation_score" = "Sequence conservation",
  "aa_diversity" = "Amino acid diversity",
  "aa_hydrophobicity" = "Hydrophobicity"
)

p2 <- ggplot(features_long, aes(x = region, y = value, fill = region)) +
  geom_boxplot(alpha = 0.7) +
  facet_wrap(~feature, scales = "free_y", labeller = labeller(feature = feature_labels)) +
  scale_fill_manual(values = c("Conserved" = "#2E86C1", "Extension" = "#E74C3C")) +
  labs(x = "", y = "Value") +
  theme_minimal(base_size = 12) +
  theme(legend.position = "none",
        strip.text = element_text(size = 12),
        axis.text.x = element_text(angle = 45, hjust = 1))

# ============================================================
# Figure 3: Relationship between species length and N-terminal extension
# ============================================================

# Compute gap ratio in the extreme N-terminal region (distance < -200) for each species
n_tail_positions <- which(distance_vec < -200)

species_extension <- apply(alignment_matrix[, n_tail_positions, drop = FALSE], 1, function(row) {
  sum(row == "-") / length(row)
})

df_species <- data.frame(
  ID = seq_names,
  Length = species_info$Length,
  N_tail_gap_ratio = species_extension
)

# Group by length
df_species <- df_species %>%
  mutate(Length_group = case_when(
    Length <= 200 ~ "Short (≤200)",
    Length <= 393 ~ "Medium (201-393)",
    Length > 393  ~ "Long (>393)"
  ))

p3 <- ggplot(df_species, aes(x = Length, y = N_tail_gap_ratio)) +
  geom_point(aes(color = Length_group), size = 2, alpha = 0.6) +
  geom_smooth(method = "lm", color = "black", se = TRUE) +
  scale_color_manual(values = c("blue", "green", "red")) +
  labs(title = "Protein length vs N-terminal extension depth",
       x = "Protein length (aa)",
       y = "Gap ratio in extreme N-terminus (distance < -200)",
       color = "Length group") +
  theme_minimal(base_size = 14) +
  theme(plot.title = element_text(hjust = 0.5))

# ============================================================
# Save all figures
# ============================================================

ggsave("Fig1_extension_profile.pdf", plot = p1, width = 10, height = 6)
ggsave("Fig2_region_features.pdf", plot = p2, width = 10, height = 6)
ggsave("Fig3_length_vs_extension.pdf", plot = p3, width = 10, height = 6)

cat("All figures generated successfully!\n")