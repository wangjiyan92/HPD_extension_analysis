================================================================
Supplementary Code for
"Structure outlasts sequence: domain emergence by terminal extension"
================================================================

1. System requirements
   - Operating system: Windows 10/11, macOS, or Linux
   - R version 4.3.2 or higher
   - Required R packages: tidyverse, Biostrings, ranger, caret, ggplot2,
     patchwork, corrplot, psych, ggpubr, gridExtra
   - No non-standard hardware is required
   - Tested on R 4.3.2 under Windows 10/11 (64-bit)

2. Installation guide
   - Install R from https://cran.r-project.org/
   - Install the required packages by running the following command in R:
     install.packages(c("tidyverse", "Biostrings", "ranger", "caret",
                        "ggplot2", "patchwork", "corrplot", "psych",
                        "ggpubr", "gridExtra"))
   - Typical install time: approximately 5-10 minutes on a normal desktop
     computer with a stable internet connection.

3. Demo
   - The script Shuffling_control_sliding_window.R contains the VOC1 and
     VOC2 nucleotide sequences embedded in the script, so it can be run
     directly without any external input file.
   - To run the demo:
     Rscript Shuffling_control_sliding_window.R
   - Expected output: a PDF file named
     "VOC1_VOC2_original_vs_shuffled.pdf" containing two sliding-window
     identity plots (original vs shuffled).
   - Expected run time: less than 1 minute on a normal desktop computer.

4. Instructions for use
   - Each R script is self-contained. Set the working directory to the
     folder containing the script and input files, then run:
     Rscript <script_name>.R
   - Input file formats and required columns are described in the header
     comments of each script.
   - Output files (CSV and PDF) are written to the working directory.

5. Script descriptions

   (1) HPD_training_data_generator.R
       - Input: avian_hpd_data.csv, avian_hpd_alignment.fasta
       - Output: hpd_training_data.csv
       - This script generates the training dataset from multi-sequence
         alignment.

   (2) HPD_extension_model_builder.R
       - Input: hpd_training_data.csv
       - Output: rf_model_extension.rds, feature_names.rds, PDF figures
       - This script trains the Random Forest model and evaluates its
         performance.

   (3) HPD_extension_visualization.R
       - Input: hpd_training_data.csv, avian_hpd_data.csv,
         avian_hpd_alignment.fasta
       - Output: 3 PDF figures (extension profile, region features,
         length vs extension)
       - This script generates the main visualizations.

   (4) Shuffling_control_sliding_window.R
       - Input: VOC1 and VOC2 nucleotide sequences (embedded)
       - Output: sliding-window identity plot (original vs shuffled)
       - This script tests whether VOC1 originated from exon shuffling.

   (5) Cross_species_validation_pisces.R
       - Input: rf_model_extension.rds, pisces_hpd_training_data.csv
       - Output: cross-species validation scatter plot (PDF)
       - This script validates the extension rule on Pisces HPD data.

   (6) Cross_protein_validation_glod4.R
       - Input: rf_model_extension.rds, glod4_training_data.csv
       - Output: cross-protein validation scatter plot (PDF)
       - This script validates the extension rule on Aves GLOD4 data.

   (7) HPD_IUPred3_aligned_plot.R
       - Input: iupred3_scores.csv (contains position, score, Length,
         VOC2_C for each isoform)
       - Output: IUPred3 disorder profile plot (PDF)
       - This script visualizes IUPred3 disorder scores of avian HPD
         isoforms aligned by VOC2_C.

   (8) multi_test_correction.R
       - Input: HPD.csv, GLOD4.csv, GLO1.csv, MCEE.csv (each containing
         Species, Length, and six trait columns)
       - Output: correlation_results_all_with_BH.csv,
         significant_results_after_BH.csv, correlation_summary_table.csv
       - This script performs Pearson correlation analyses across
         4 proteins x 6 traits and applies Benjamini-Hochberg correction
         for multiple testing.

6. Application to other species/proteins
   The same analysis logic was applied to other datasets (Insecta,
   Mammalia, Plantae, Mammalian GLOD4, GLO1, MCEE) by replacing the input
   file names and reference IDs. The scripts above demonstrate the
   complete analytical workflow.

7. Reproduction instructions
   To reproduce the main quantitative results in the manuscript:
   1. Run HPD_training_data_generator.R to generate hpd_training_data.csv.
   2. Run HPD_extension_model_builder.R to train the Random Forest model
      and generate the performance metrics (R2 and RMSE).
   3. Run HPD_extension_visualization.R to reproduce the main figures.
   4. Run Shuffling_control_sliding_window.R to reproduce the
      sliding-window control analysis.
   5. Run Cross_species_validation_pisces.R and
      Cross_protein_validation_glod4.R to reproduce the cross-species and
      cross-protein validation results.
   6. Run HPD_IUPred3_aligned_plot.R to reproduce the IUPred3 disorder
      profile.
   7. Run multi_test_correction.R to reproduce the multiple testing
      correction results.

8. Contact
   For questions regarding the code, please contact:
   wangjiyan@ustc.edu.cn
