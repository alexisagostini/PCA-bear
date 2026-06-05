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
library(ggplot2)

# Chargement PCA
eigenvec <- read.table("/home/alexis/data/Ours/results/ours_unique_pca.eigenvec",
                       header=TRUE, comment.char="")
colnames(eigenvec)[1] <- "IID"

eigenval <- scan("/home/alexis/data/Ours/results/ours_unique_pca.eigenval")
pct <- round(eigenval / sum(eigenval) * 100, 2)

# Chargement métadonnées — jointure directe sur run_accession
meta <- read.csv("/home/alexis/data/Ours/results/SRA_Runtable_final.csv",
                 sep=";", header=TRUE, stringsAsFactors=FALSE)

meta_sub <- meta[, c("run_accession",
                      "geographic.location..country.and.or.sea.",
                      "Continent.location",
                      "sex",
                      "collection.date",
                      "submitter.id")]
colnames(meta_sub) <- c("IID", "country", "continent", "sex", "year", "sample_id")

df <- merge(eigenvec, meta_sub, by="IID", all.x=TRUE)

# Vérification jointure
cat("Lignes avec country NA :", sum(is.na(df$country)), "\n")
cat("Pays présents :", paste(unique(df$country), collapse=", "), "\n")

# Plot PC1 vs PC2 coloré par pays
p1 <- ggplot(df, aes(x=PC1, y=PC2, color=country, label=sample_id)) +
  geom_point(size=3, alpha=0.85) +
  geom_text(size=2.2, vjust=-0.8, hjust=0.5, show.legend=FALSE) +
  labs(
    title="PCA Ursus arctos — 57 individus uniques",
    x=paste0("PC1 (", pct[1], "%)"),
    y=paste0("PC2 (", pct[2], "%)"),
    color="Pays"
  ) +
  theme_bw(base_size=13) +
  theme(legend.position="right")

# Plot PC1 vs PC3
p2 <- ggplot(df, aes(x=PC1, y=PC3, color=country, label=sample_id)) +
  geom_point(size=3, alpha=0.85) +
  geom_text(size=2.2, vjust=-0.8, hjust=0.5, show.legend=FALSE) +
  labs(
    title="PCA Ursus arctos — PC1 vs PC3",
    x=paste0("PC1 (", pct[1], "%)"),
    y=paste0("PC3 (", pct[3], "%)"),
    color="Pays"
  ) +
  theme_bw(base_size=13) +
  theme(legend.position="right")

ggsave("/home/alexis/data/Ours/results/pca/pca_PC1_PC2.pdf", p1, width=12, height=8)
ggsave("/home/alexis/data/Ours/results/pca/pca_PC1_PC3.pdf", p2, width=12, height=8)
cat("Plots sauvegardés dans /home/alexis/data/Ours/results/pca/\n")
```
