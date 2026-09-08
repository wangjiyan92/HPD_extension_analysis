# Load required packages
library(tidyverse)
library(ggplot2)

# Read data
data <- read.csv("iupred3_scores.csv", stringsAsFactors = FALSE)

# Calculate distance relative to VOC2_C
data <- data %>%
  mutate(distance = position - VOC2_C)

# Custom colors (adjustable)
my_colors <- c(
  "140" = "#1f77b4",   # dark blue
  "180" = "#ff7f0e",   # orange
  "187" = "#2ca02c",   # green
  "312" = "#d62728",   # red
  "354" = "#9467bd",   # purple
  "363" = "#8c564b",   # brown
  "393" = "#e377c2"    # pink
)

# Plot
p <- ggplot(data, aes(x = distance, y = score, color = as.factor(Length), group = Length)) +
  geom_line(size = 1.2) +
  geom_hline(yintercept = 0.5, linetype = "dashed", color = "gray50", alpha = 0.7) +
  geom_vline(xintercept = 0, linetype = "dashed", color = "red", alpha = 0.5) +
  annotate("text", x = 0, y = 0.95, label = "VOC2_C", color = "red", hjust = 0.5, size = 4) +
  annotate("text", x = -50, y = 0.95, label = "0.5 threshold", color = "gray50", hjust = 0, size = 4) +
  ylim(0, 1) +
  labs(
    title = "IUPred3 disorder scores aligned by VOC2_C",
    x = "Distance from VOC2_C (negative = N-terminus, positive = C-terminus)",
    y = "Disorder score",
    color = "Length (aa)"
  ) +
  theme_minimal(base_size = 14) +
  theme(
    legend.position = "right",
    plot.title = element_text(hjust = 0.5)
  ) +
  scale_color_manual(values = my_colors)

print(p)

# Save as PDF
ggsave("iupred3_aligned_by_VOC2_C.pdf", plot = p, width = 10, height = 6)