# PCA — Ursus arctos

## preparation and pca files

```
#!/bin/bash
#SBATCH --job-name=ours_pca
#SBATCH --output=/home/alexis/data/Ours/results/pca/final-pca/ours_pca_%j.log
#SBATCH --error=/home/alexis/data/Ours/results/pca/final-pca/ours_pca_%j.err

# Activation de l'environnement conda
source activate ours_pca

# Déduplication
vcftools --gzvcf /home/alexis/data/Ours/results/final_variants.clean.vcf.gz \
  --keep /home/alexis/data/Ours/results/samples_to_keep.txt \
  --recode --stdout | bgzip -c > /home/alexis/data/Ours/results/final_variants.clean.dedup.vcf.gz
tabix -p vcf /home/alexis/data/Ours/results/final_variants.clean.dedup.vcf.gz

# Vérification
bcftools query -l /home/alexis/data/Ours/results/final_variants.clean.dedup.vcf.gz | wc -l

mkdir -p /home/alexis/data/Ours/results/pca/final-pca
cd /home/alexis/data/Ours/results/pca/final-pca

# Filtrage MAF 0.05
vcftools --gzvcf /home/alexis/data/Ours/results/final_variants.clean.dedup.vcf.gz \
  --maf 0.05 \
  --recode --stdout | bgzip -c > output.maf05.vcf.gz
tabix -p vcf output.maf05.vcf.gz

# VCF → BED
plink --vcf output.maf05.vcf.gz \
  --allow-extra-chr \
  --make-bed \
  --out data

# Pruning LD
plink --bfile data \
  --allow-extra-chr \
  --indep-pairwise 50 5 0.2 \
  --out prune

# Appliquer le pruning
plink --bfile data \
  --allow-extra-chr \
  --extract prune.prune.in \
  --make-bed \
  --out data_pruned

# PCA
plink --bfile data_pruned \
  --allow-extra-chr \
  --pca \
  --out PCA

# ROH
plink --bfile data_pruned \
  --allow-extra-chr \
  --homozyg \
  --homozyg-snp 50 \
  --homozyg-kb 1000 \
  --homozyg-density 50 \
  --homozyg-gap 1000 \
  --out ROH
```
## PCA
```R
#!/usr/bin/env Rscript
# ============================================================
# Génération des plots PCA et ROH à partir du VCF
# à partir de : /home/alexis/data/Ours/results/final_variants.clean.vcf.gz
# ============================================================

library(ggplot2)
library(dplyr)
library(tidyr)

# etup paths
VCF_PATH        <- "/home/alexis/data/Ours/results/final_variants.clean.vcf.gz"
METADATA_PATH   <- "/home/alexis/data/Ours/results/SRA_Runtable_final.csv"
OUTDIR          <- "/home/alexis/data/Ours/results/pca/final-pca/Png"
PLINK_DIR       <- "/home/alexis/data/Ours/results/pca/final-pca"

dir.create(OUTDIR, showWarnings = FALSE, recursive = TRUE)

# 1. Load metadata 
cat("Loading metadata...\n")
meta <- read.csv(METADATA_PATH, sep = ";", stringsAsFactors = FALSE, check.names = FALSE)

# Extract relevant columns
meta_clean <- meta %>%
  select(
    run_accession,
    `Continent location`,
    `geographic location (country and/or sea)`
  ) %>%
  rename(
    Sample = run_accession,
    Continent = `Continent location`,
    Country = `geographic location (country and/or sea)`
  ) %>%
  filter(!is.na(Sample)) %>%
  distinct()

cat("Metadata samples:", nrow(meta_clean), "\n")
cat("Continents:", paste(unique(meta_clean$Continent), collapse = ", "), "\n")

# 2. Load PCA results from PLINK 
cat("\nLoading PCA results...\n")

eigenvec_path <- file.path(PLINK_DIR, "PCA.eigenvec")
eigenval_path <- file.path(PLINK_DIR, "PCA.eigenval")

if (!file.exists(eigenvec_path) || !file.exists(eigenval_path)) {
  stop(sprintf("PCA files not found:\n  %s\n  %s", eigenvec_path, eigenval_path))
}

eigenvec <- read.table(eigenvec_path, header = FALSE, stringsAsFactors = FALSE)
eigenval <- scan(eigenval_path)

# PLINK eigenvec format: FID IID PC1 PC2 PC3 ... (no header)
# Assign column names: first two columns are FID, IID, then PC1-PC18
n_pcs <- ncol(eigenvec) - 2
colnames(eigenvec) <- c("FID", "IID", paste0("PC", 1:n_pcs))

pca_df <- eigenvec %>%
  select(FID, IID, PC1, PC2, PC3) %>%
  rename(Sample = IID)

# Calculate variance explained
total_var <- sum(eigenval)
var_pct <- round(100 * eigenval / total_var, 2)

cat("PC1:", var_pct[1], "%\n")
cat("PC2:", var_pct[2], "%\n")
cat("PC3:", var_pct[3], "%\n")

# 3. Merge PCA + metadata 
cat("\nMerging PCA and metadata...\n")
pca_annot <- pca_df %>%
  left_join(meta_clean, by = "Sample")

# Check for unmatched samples
n_unmatched <- sum(is.na(pca_annot$Continent))
cat("Unmatched samples (no metadata):", n_unmatched, "/", nrow(pca_annot), "\n")

# Fill missing values
pca_annot$Continent[is.na(pca_annot$Continent)] <- "Unknown"
pca_annot$Country[is.na(pca_annot$Country)]     <- "Unknown"

# 4. Plot PCA (PC1 vs PC2) colored by Continent 
cat("\nPlotting PCA (PC1 vs PC2, colored by Continent)...\n")

p_pca_continent <- ggplot(pca_annot, aes(x = PC1, y = PC2, color = Continent)) +
  geom_point(size = 3, alpha = 0.7) +
  labs(
    title = "PCA — PC1 vs PC2",
    x = sprintf("PC1 (%.2f%%)", var_pct[1]),
    y = sprintf("PC2 (%.2f%%)", var_pct[2]),
    color = "Continent"
  ) +
  theme_bw(base_size = 12) +
  theme(
    legend.position = "right",
    plot.title = element_text(hjust = 0.5, face = "bold")
  )

ggsave(
  file.path(OUTDIR, "PCA_PC1_PC2_continent.png"),
  p_pca_continent,
  width = 10,
  height = 8,
  dpi = 150
)
cat("✓ Saved:", file.path(OUTDIR, "PCA_PC1_PC2_continent.png"), "\n")

# 5. Load ROH results from PLINK
cat("\nLoading ROH results...\n")

roh_path <- file.path(PLINK_DIR, "ROH.hom.indiv")

if (!file.exists(roh_path)) {
  stop(sprintf("ROH file not found: %s", roh_path))
}

roh_df <- read.table(roh_path, header = TRUE, stringsAsFactors = FALSE)

# Rename column to match metadata
colnames(roh_df)[colnames(roh_df) == "IID"] <- "Sample"

# Calculate FROH (Fraction of Runs of Homozygosity)
# Genome size estimation (use actual genome size if known)
# For Ursus arctos, approximate genome size ~2,400 Mb
GENOME_SIZE_KB <- 2400000  # 2.4 Gb in kb

roh_df <- roh_df %>%
  mutate(
    FROH = KB / GENOME_SIZE_KB  # Fraction of genome in ROH
  )

cat("ROH summary:\n")
cat("  Mean FROH:", round(mean(roh_df$FROH, na.rm = TRUE), 4), "\n")
cat("  Mean NSEG:", round(mean(roh_df$NSEG, na.rm = TRUE), 2), "\n")

# Merge ROH + metadata
roh_annot <- roh_df %>%
  left_join(meta_clean, by = "Sample")

roh_annot$Continent[is.na(roh_annot$Continent)] <- "Unknown"

# 6. Plot ROH (FROH by sample, colored by Continent) 
cat("\nPlotting ROH (FROH, colored by Continent)...\n")

# Sort by FROH for better visualization
roh_annot_sorted <- roh_annot %>%
  arrange(FROH)

p_roh_froh <- ggplot(roh_annot_sorted, aes(x = reorder(Sample, FROH), y = FROH, fill = Continent)) +
  geom_bar(stat = "identity", alpha = 0.8) +
  labs(
    title = "Runs of Homozygosity (ROH) — FROH by Sample",
    x = "Sample",
    y = "FROH (Fraction of Genome in ROH)",
    fill = "Continent"
  ) +
  theme_bw(base_size = 10) +
  theme(
    axis.text.x = element_text(angle = 90, vjust = 0.5, hjust = 1, size = 8),
    plot.title = element_text(hjust = 0.5, face = "bold"),
    legend.position = "right"
  )

ggsave(
  file.path(OUTDIR, "ROH_FROH.png"),
  p_roh_froh,
  width = 16,
  height = 8,
  dpi = 150
)
```
