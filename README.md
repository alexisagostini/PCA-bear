# PCA — Ursus arctos

## 1. Missingness

```bash
plink2 \
  --vcf /home/alexis/data/Ours/ours_filtered.vcf.gz \
  --missing \
  --out /home/alexis/data/Ours/results/ours_missing \
  --allow-extra-chr
```

## 2. Sample filtering

```bash
echo "ERR16661338
SRR7408244" > /home/alexis/data/Ours/results/samples_to_remove.txt

plink2 \
  --vcf /home/alexis/data/Ours/ours_filtered.vcf.gz \
  --remove /home/alexis/data/Ours/results/samples_to_remove.txt \
  --make-pgen \
  --out /home/alexis/data/Ours/results/ours_filtered_clean \
  --allow-extra-chr
```

> 295 samples retained, 166,465 variants

## 3. PCA (SLURM)

```bash
#!/bin/bash
#SBATCH --job-name=pca_ours
#SBATCH --mem=50G
#SBATCH --cpus-per-task=2
#SBATCH --time=02:00:00
#SBATCH --output=/home/alexis/data/Ours/results/pca_%j.log

plink2 \
  --pfile /home/alexis/data/Ours/results/ours_filtered_clean \
  --pca approx 10 \
  --out /home/alexis/data/Ours/results/ours_pca \
  --allow-extra-chr
```

> Output: `ours_pca.eigenvec` / `ours_pca.eigenval`

## 4. Visualization in R

```r
library(tidyverse)
library(ggrepel)

eigenvec <- read_table("/home/alexis/data/Ours/results/ours_pca.eigenvec")
eigenval <- read_lines("/home/alexis/data/Ours/results/ours_pca.eigenval") |> as.numeric()
pct_var  <- round(eigenval / sum(eigenval) * 100, 2)

meta <- read_delim("/home/alexis/data/Ours/SRA_Runtable_final.csv", delim = ";") |>
  rename(
    country   = `geographic location (country and/or sea)`,
    continent = `Continent location`
  ) |>
  select(run_accession, country, continent)

pca_meta <- eigenvec |>
  left_join(meta, by = c("#IID" = "run_accession"))

ggplot(pca_meta, aes(x = PC1, y = PC2, color = country, shape = continent)) +
  geom_point(size = 2.5, alpha = 0.8) +
  geom_text_repel(aes(label = country), size = 3, max.overlaps = 20,
                  segment.size = 0.3, segment.alpha = 0.5) +
  scale_color_brewer(palette = "Set2") +
  labs(
    title   = "PCA — Ursus arctos",
    x       = paste0("PC1 (", pct_var[1], "%)"),
    y       = paste0("PC2 (", pct_var[2], "%)"),
    color   = "Country",
    shape   = "Continent",
    caption = paste0("n = ", nrow(pca_meta), " samples | 166,465 variants")
  ) +
  theme_bw()

ggsave("/home/alexis/data/Ours/results/pca_plot_labeled.png", width = 14, height = 10, dpi = 300)
```

> Note: 34 Irish samples (SRR7408xxx) are absent — filtered out by PLINK due to empty genotype data.
