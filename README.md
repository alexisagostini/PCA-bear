# PCA — Ursus arctos

## preparation and pca

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
