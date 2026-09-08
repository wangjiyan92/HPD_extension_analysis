================================================================
Supplementary Code for
"Structure outlasts sequence: domain emergence by terminal extension"
================================================================

1. Requirements
   - R version 4.3.2 or higher
   - Required R packages: tidyverse, Biostrings, ranger, caret, ggplot2, patchwork, corrplot, psych, ggpubr, gridExtra
   - Install with: install.packages(c("tidyverse", "Biostrings", "ranger", "caret", "ggplot2", "patchwork", "corrplot", "psych", "ggpubr", "gridExtra"))

2. Script descriptions

   (1) HPD_training_data_generator.R
       - Input: avian_hpd_data.csv, avian_hpd_alignment.fasta
       - Output: hpd_training_data.csv
       - This script generates the training dataset from multi-sequence alignment.

   (2) HPD_extension_model_builder.R
       - Input: hpd_training_data.csv
       - Output: rf_model_extension.rds, feature_names.rds, PDF figures
       - This script trains the Random Forest model and evaluates its performance.

   (3) HPD_extension_visualization.R
       - Input: hpd_training_data.csv, avian_hpd_data.csv, avian_hpd_alignment.fasta
       - Output: 3 PDF figures (extension profile, region features, length vs extension)
       - This script generates the main visualizations.

   (4) Shuffling_control_sliding_window.R
       - Input: VOC1 and VOC2 nucleotide sequences (embedded)
       - Output: Sliding-window identity plot (original vs shuffled)
       - This script tests whether VOC1 originated from exon shuffling.

   (5) Cross_species_validation_pisces.R
       - Input: rf_model_extension.rds, pisces_hpd_training_data.csv
       - Output: Cross-species validation scatter plot (PDF)
       - This script validates the extension rule on Pisces HPD data.

   (6) Cross_protein_validation_glod4.R
       - Input: rf_model_extension.rds, glod4_training_data.csv
       - Output: Cross-protein validation scatter plot (PDF)
       - This script validates the extension rule on Aves GLOD4 data.

   (7) HPD_IUPred3_aligned_plot.R
       - Input: iupred3_scores.csv (contains position, score, Length, VOC2_C for each isoform)
       - Output: IUPred3 disorder profile plot (PDF)
       - This script visualizes IUPred3 disorder scores of avian HPD isoforms aligned by VOC2_C.

   (8) multi_test_correction.R
       - Input: HPD.csv, GLOD4.csv, GLO1.csv, MCEE.csv (each containing Species, Length, and six trait columns)
       - Output: correlation_results_all_with_BH.csv, significant_results_after_BH.csv, correlation_summary_table.csv
       - This script performs Pearson correlation analyses across 4 proteins × 6 traits and applies Benjamini-Hochberg correction for multiple testing.

3. Application to other species/proteins
   The same analysis logic was applied to other datasets (Insecta, Mammalia,
   Plantae, Mammalian GLOD4, GLO1, MCEE) by replacing the input file names and reference IDs.
   The scripts above demonstrate the complete analytical workflow.

4. Contact
   For questions regarding the code, please contact: wangjiyan@ustc.edu.cn