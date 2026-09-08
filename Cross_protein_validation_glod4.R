# ============================================================
# Cross-protein validation: Aves-trained HPD model applied to Aves GLOD4 data
# Purpose: Test whether the N-terminal extension rule observed in HPD
#          applies to GLOD4 within the same protein family
# Input:
#   - rf_model_extension.rds (Aves-trained HPD random forest model)
#   - glod4_training_data.csv (Aves GLOD4 training data, same format as HPD)
# Output:
#   - Predicted vs actual scatter plot (PDF)
#   - Evaluation metrics (R², RMSE)
# ============================================================

library(tidyverse)
library(ranger)
library(ggplot2)
library(caret)

# ============================================================
# 1. Load Aves-trained HPD model
# ============================================================

model <- readRDS("rf_model_extension.rds")
cat("Aves HPD model loaded\n")

# Load feature names used during training
feature_names <- readRDS("feature_names.rds")
cat("Training features:", paste(feature_names, collapse = ", "), "\n")

# ============================================================
# 2. Load Aves GLOD4 data
# ============================================================

glod4_data <- read.csv("glod4_training_data.csv", stringsAsFactors = FALSE)

# Check data
cat("Aves GLOD4 data rows:", nrow(glod4_data), "\n")

# Ensure categorical variables are factors
factor_cols <- c("aa_group", "aa_charge")
for (col in factor_cols) {
  if (col %in% colnames(glod4_data)) {
    glod4_data[[col]] <- as.factor(glod4_data[[col]])
  }
}

# ============================================================
# 3. Feature preparation (identical to training set)
# ============================================================

# Extract features
X_glod4 <- glod4_data[, feature_names, drop = FALSE]

# Handle missing values (median imputation for numeric, consistent with training)
for (col in names(X_glod4)) {
  if (is.numeric(X_glod4[[col]]) && any(is.na(X_glod4[[col]]))) {
    med_val <- median(X_glod4[[col]], na.rm = TRUE)
    X_glod4[[col]][is.na(X_glod4[[col]])] <- med_val
    cat("Imputed missing values in", col, "with median:", med_val, "\n")
  }
}

# Extract actual values
y_true <- glod4_data$extension_potential

# ============================================================
# 4. Predict GLOD4 data using Aves-trained HPD model
# ============================================================

y_pred <- predict(model, data = X_glod4)$predictions

# ============================================================
# 5. Calculate evaluation metrics
# ============================================================

rmse <- RMSE(y_pred, y_true)
r2 <- R2(y_pred, y_true)

cat("\n=== Cross-protein prediction evaluation ===\n")
cat("Model: Aves HPD -> Target: Aves GLOD4\n")
cat("R² =", round(r2, 4), "\n")
cat("RMSE =", round(rmse, 4), "\n")

# ============================================================
# 6. Generate scatter plot (with diagonal line)
# ============================================================

plot_data <- data.frame(
  Actual = y_true,
  Predicted = y_pred
)

p <- ggplot(plot_data, aes(x = Actual, y = Predicted)) +
  geom_point(alpha = 0.3, color = "darkorange", size = 1.2) +
  geom_abline(slope = 1, intercept = 0, linetype = "dashed", color = "black", size = 0.8) +
  labs(
    title = "Cross-protein validation: Aves-trained HPD model applied to Aves GLOD4",
    x = "Actual extension potential (Aves GLOD4)",
    y = "Predicted extension potential (by Aves HPD model)"
  ) +
  annotate("text", x = 0.1, y = 0.85, 
           label = paste0("R² = ", round(r2, 3), "\nRMSE = ", round(rmse, 3)),
           hjust = 0, size = 4) +
  coord_fixed(ratio = 1) +
  theme_minimal(base_size = 14) +
  theme(
    plot.title = element_text(hjust = 0.5),
    panel.grid.minor = element_blank()
  )

# Save as PDF
ggsave("cross_protein_validation_glod4.pdf", plot = p, width = 8, height = 7)

# Also save as PNG
ggsave("cross_protein_validation_glod4.png", plot = p, width = 8, height = 7, dpi = 300)

cat("\nScatter plot saved: cross_protein_validation_glod4.pdf\n")