# ============================================================
# Destructive validation: shuffled VOC1 vs VOC2 sliding-window analysis
# Purpose: Test whether the observed high identity between VOC1 and VOC2
#          is caused by sequence composition bias
# Method: Randomly shuffle the VOC1 nucleotide sequence (preserving base composition),
#         then repeat the same analysis against VOC2, comparing the original
#         and shuffled curves
# Output: PDF containing two subplots (original vs shuffled)
# ============================================================

# Load required packages
library(Biostrings)
library(ggplot2)
library(patchwork)  # For combining plots

# ============================================================
# 1. Define sequences (VOC1 and VOC2 provided by user)
# ============================================================

voc1_seq <- "CACTTCCACTCTGTGACCTTCTGGGTTGGCAACGCCAAGCAGGCCGCGTCATTCTACTGCAGCAAGATGGGCTTTGAACCTCTAGCCTACAGGGGCCTGGAGACCGGTTCCCGGGAGGTGGTCAGCCATGTAATCAAACAAGGGAAGATTGTGTTTGTCCTCTCCTCAGCGCTCAACCCCTGGAACAAAGAGATGGGCGATCACCTGGTGAAACACGGTGACGGAGTGAAGGACATTGCGTTCGAGGTGGAAGATTGTGACTACATCGTGCAGAAAGCACGGGAACGGGGCGCCAAAATCATGCGGGAGCCCTGGGTAGAGCAAGACAAGTTTGGGAAGGTGAAGTTTGCTGTGCTGCAGACGTATGGGGACACCACACACACCCTGGTGGAGAAG"

voc2_seq <- "ATGATCGACCACATTGTGGGAAACCAGCCTGATCAGGAGATGGTGTCCGCCTCCGAATGGTACCTGAAAAACCTGCAGTTCCACCGCTTCTGGTCCGTGGATGACACGCAGGTGCACACGGAATATAGCTCTCTGCGATCCATTGTGGTGGCCAACTATGAAGAGTCCATCAAGATGCCCATCAATGAGCCAGCGCCTGGCAAGAAGAAGTCCCAGATCCAGGAATATGTGGACTATAACGGGGGCGCTGGGGTCCAGCACATCGCTCTCAAGACCGAAGACATCATCACAGCGATTCGCCACTTGAGAGAGAGAGGCCTGGAGTTCTTATCTGTTCCCTCCACGTACTACAAACAACTGCGGGAGAAGCTGAAGACGGCCAAGATCAAGGTGAAGGAGAACATTGATGCCCTGGAGGAGCTGAAAATCCTGGTGGACTACGACGAGAAAGGCTACCTCCTGCAGATCTTCACCAAA"

# ============================================================
# 2. Define sliding-window analysis function
# ============================================================

sliding_window_analysis <- function(seq1, seq2, window_size = 27, step = 3) {
  # Convert sequences to DNAString
  dna1 <- DNAString(seq1)
  dna2 <- DNAString(seq2)
  
  # Global alignment
  aln <- pairwiseAlignment(dna1, dna2, type = "global",
                           gapOpening = 10, gapExtension = 4,
                           substitutionMatrix = "BLOSUM62")
  aln1 <- as.character(pattern(aln))
  aln2 <- as.character(subject(aln))
  
  total_len <- nchar(aln1)
  if (total_len < window_size) {
    stop("Aligned sequence length is shorter than window size")
  }
  
  positions <- seq(1, total_len - window_size + 1, by = step)
  similarities <- numeric(length(positions))
  
  for (i in seq_along(positions)) {
    start <- positions[i]
    end <- start + window_size - 1
    win1 <- substr(aln1, start, end)
    win2 <- substr(aln2, start, end)
    chars1 <- strsplit(win1, "")[[1]]
    chars2 <- strsplit(win2, "")[[1]]
    valid <- chars1 != "-" & chars2 != "-"
    if (sum(valid) > 0) {
      matches <- sum(chars1[valid] == chars2[valid])
      similarities[i] <- matches / sum(valid)
    } else {
      similarities[i] <- NA
    }
  }
  
  data.frame(Position = positions, Identity = similarities)
}

# ============================================================
# 3. Original VOC1 vs VOC2 analysis
# ============================================================

df_original <- sliding_window_analysis(voc1_seq, voc2_seq)
df_original_clean <- df_original[!is.na(df_original$Identity), ]

# ============================================================
# 4. Destructive validation: shuffled VOC1 vs VOC2 (10 permutations averaged)
# ============================================================

set.seed(12345)  # Fixed random seed for reproducibility

n_permutations <- 10  # 10 permutations to reduce random fluctuation

# Store all shuffled results
all_shuffled <- list()

for (p in 1:n_permutations) {
  # Shuffle VOC1 sequence (preserving base composition)
  chars <- strsplit(voc1_seq, "")[[1]]
  shuffled_chars <- sample(chars, length(chars), replace = FALSE)
  shuffled_seq <- paste(shuffled_chars, collapse = "")
  
  # Run sliding-window analysis
  df_shuffled <- sliding_window_analysis(shuffled_seq, voc2_seq)
  df_shuffled_clean <- df_shuffled[!is.na(df_shuffled$Identity), ]
  
  all_shuffled[[p]] <- df_shuffled_clean
}

# Combine all shuffled results and average by position
positions_common <- Reduce(intersect, lapply(all_shuffled, function(x) x$Position))
if (length(positions_common) == 0) {
  positions_common <- all_shuffled[[1]]$Position
}

# Extract Identity values at common positions and compute mean and SD
identity_matrix <- matrix(NA, nrow = length(positions_common), ncol = n_permutations)
for (p in 1:n_permutations) {
  idx <- match(positions_common, all_shuffled[[p]]$Position)
  identity_matrix[, p] <- all_shuffled[[p]]$Identity[idx]
}
mean_identity <- rowMeans(identity_matrix, na.rm = TRUE)
sd_identity <- apply(identity_matrix, 1, sd, na.rm = TRUE)

df_shuffled_summary <- data.frame(
  Position = positions_common,
  Mean_Identity = mean_identity,
  SD_Identity = sd_identity
)

# ============================================================
# 5. Plot comparison
# ============================================================

# Plot 1: Original curve
p1 <- ggplot(df_original_clean, aes(x = Position, y = Identity)) +
  geom_line(color = "steelblue", size = 1) +
  geom_hline(yintercept = 0.3, linetype = "dashed", color = "red", alpha = 0.6) +
  annotate("text", x = min(df_original_clean$Position), y = 0.32,
           label = "30% threshold", hjust = 0, size = 3, color = "red") +
  labs(title = "Original: VOC1 vs VOC2",
       x = "Alignment position (bp)",
       y = "Local nucleotide identity") +
  ylim(0, 1) +
  theme_minimal(base_size = 12) +
  theme(plot.title = element_text(hjust = 0.5))

# Plot 2: Shuffled curve (mean ± SD)
p2 <- ggplot(df_shuffled_summary, aes(x = Position, y = Mean_Identity)) +
  geom_ribbon(aes(ymin = Mean_Identity - SD_Identity, ymax = Mean_Identity + SD_Identity),
              fill = "gray80", alpha = 0.5) +
  geom_line(color = "darkorange", size = 1) +
  geom_hline(yintercept = 0.3, linetype = "dashed", color = "red", alpha = 0.6) +
  annotate("text", x = min(df_shuffled_summary$Position), y = 0.32,
           label = "30% threshold", hjust = 0, size = 3, color = "red") +
  labs(title = "Shuffled VOC1 vs VOC2 (mean ± SD, 10 permutations)",
       x = "Alignment position (bp)",
       y = "Local nucleotide identity") +
  ylim(0, 1) +
  theme_minimal(base_size = 12) +
  theme(plot.title = element_text(hjust = 0.5))

# Combine two plots (side by side)
p_combined <- p1 + p2 + plot_annotation(
  title = "Destructive validation: VOC1 shuffled control experiment",
  subtitle = "If the shuffled curve remains above 30%, the high identity reflects compositional bias"
)

# ============================================================
# 6. Save results
# ============================================================

ggsave("VOC1_VOC2_original_vs_shuffled.pdf", plot = p_combined, width = 14, height = 6)
cat("Comparison plot saved: VOC1_VOC2_original_vs_shuffled.pdf\n")

# Also save individual plots
ggsave("VOC1_VOC2_original.pdf", plot = p1, width = 7, height = 5)
ggsave("VOC1_VOC2_shuffled.pdf", plot = p2, width = 7, height = 5)

# ============================================================
# 7. Statistical summary output
# ============================================================

cat("\n=== Statistical summary ===\n")
cat("Original analysis: valid windows =", nrow(df_original_clean), "\n")
cat("Original analysis: mean identity =", round(mean(df_original_clean$Identity), 3), "\n")
cat("Original analysis: max identity =", round(max(df_original_clean$Identity), 3), "\n")
cat("Original analysis: min identity =", round(min(df_original_clean$Identity), 3), "\n\n")

cat("Shuffled analysis (10 permutations averaged):\n")
cat("Valid windows =", nrow(df_shuffled_summary), "\n")
cat("Mean identity =", round(mean(df_shuffled_summary$Mean_Identity), 3), "\n")
cat("Maximum mean identity =", round(max(df_shuffled_summary$Mean_Identity), 3), "\n")
cat("Minimum mean identity =", round(min(df_shuffled_summary$Mean_Identity), 3), "\n")