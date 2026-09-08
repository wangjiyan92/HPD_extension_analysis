# ============================================================
# Cross-species validation: Aves-trained HPD model applied to Pisces HPD data
# Purpose: Test the cross-species generality of the N-terminal extension rule
# Input:
#   - rf_model_extension.rds (Aves-trained random forest model)
#   - pisces_hpd_training_data.csv (Pisces training data, same format as Aves)
# Output:
#   - Predicted vs actual scatter plot (PDF)
#   - Evaluation metrics (R², RMSE)
# ============================================================

library(tidyverse)
library(ranger)
library(ggplot2)
library(caret)

# ============================================================
# 1. Load Aves-trained model
# ============================================================

model <- readRDS("rf_model_extension.rds")
cat("Aves model loaded\n")

# Load feature names used during training
feature_names <- readRDS("feature_names.rds")
cat("Training features:", paste(feature_names, collapse = ", "), "\n")

# ============================================================
# 2. Load Pisces HPD data
# ============================================================

pisces_data <- read.csv("pisces_hpd_training_data.csv", stringsAsFactors = FALSE)

# Check data
cat("Pisces data rows:", nrow(pisces_data), "\n")

# Ensure categorical variables are factors
factor_cols <- c("aa_group", "aa_charge")
for (col in factor_cols) {
  if (col %in% colnames(pisces_data)) {
    pisces_data[[col]] <- as.factor(pisces_data[[col]])
  }
}

# ============================================================
# 3. Feature preparation (identical to training set)
# ============================================================

# Extract features
X_pisces <- pisces_data[, feature_names, drop = FALSE]

# Handle missing values (median imputation for numeric, consistent with training)
for (col in names(X_pisces)) {
  if (is.numeric(X_pisces[[col]]) && any(is.na(X_pisces[[col]]))) {
    med_val <- median(X_pisces[[col]], na.rm = TRUE)
    X_pisces[[col]][is.na(X_pisces[[col]])] <- med_val
    cat("Imputed missing values in", col, "with median:", med_val, "\n")
  }
}

# Extract actual values
y_true <- pisces_data$extension_potential

# ============================================================
# 4. Predict Pisces HPD data using Aves-trained model
# ============================================================

y_pred <- predict(model, data = X_pisces)$predictions

# ============================================================
# 5. Calculate evaluation metrics
# ============================================================

rmse <- RMSE(y_pred, y_true)
r2 <- R2(y_pred, y_true)

cat("\n=== Cross-species prediction evaluation ===\n")
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
    title = "Cross-species validation: Aves-trained HPD model applied to Pisces HPD",
    x = "Actual extension potential (Pisces)",
    y = "Predicted extension potential (by Aves model)"
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
ggsave("cross_species_validation_pisces.pdf", plot = p, width = 8, height = 7)

# Also save as PNG
ggsave("cross_species_validation_pisces.png", plot = p, width = 8, height = 7, dpi = 300)

cat("\nScatter plot saved: cross_species_validation_pisces.pdf\n")