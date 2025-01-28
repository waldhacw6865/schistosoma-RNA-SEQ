library(tidyverse)
library(DESeq2)

counts_df <- read_tsv("/data/users/waldhacw6865/Sm_RNA-seq_Project/Notebook/counts/counts.tsv",
                      comment = "#") |>
             mutate(across(where(is.numeric), as.integer))

counts_summary <- counts_df |>
    select(Geneid, contains('.bam')) |>
    rename_with(~str_remove(., "/dedup/star.bam:"), everything()) |>
    rowwise() |>
    mutate(total_counts = sum(c_across(where(is.numeric)), na.rm = T)) |>
    filter(total_counts >= 10)

sample_summary <- counts_df |>
    select(Geneid, contains('.bam')) |>
    rename_with(~str_remove(., "/dedup/star.bam:"), everything()) |>
    pivot_longer(-Geneid, names_to = 'sample', values_to = 'count') |>
    filter(count > 0) |>
    group_by(Geneid) |>
    tally() |>
    filter(n <= 3)

genes_to_remove = sample_summary$Geneid

counts_filt <- counts_summary |>
    filter(!Geneid %in% genes_to_remove) |>
    arrange(Geneid) |>
    select(-total_counts)
counts_filt

counts_m <- counts_filt |>
    select(-Geneid) |>
    as.matrix()
rownames(counts_m) <- counts_filt$Geneid
counts_m

dists <- dist(t(counts_m))
dists

dists_df <- as.matrix(dists) |>
    as_tibble(rownames = 'sample')

dist_plot <- dists_df |>
    pivot_longer(-sample, names_to = 'comp', values_to = 'dist') |>
    ggplot(aes(x = sample, y = comp, fill = dist)) +
    geom_tile() +
    scale_fill_viridis_c() +
    coord_equal() +
    NULL
dist_plot

pca_fit <- t(log10(counts_m + 1)) |> 
  prcomp(scale = TRUE)
pca_fit

library(broom)

metadata <- data.frame(sample_id = colnames(counts_m)) |>
    mutate(tissue = str_sub(sample_id, 12, 14),
           rep = str_sub(sample_id, 17, 17))
rownames(metadata) <- metadata$sample_id
metadata <- select(metadata, -sample_id)
metadata

all(rownames(metadata) == colnames(counts_m))

dds <- DESeqDataSetFromMatrix(countData = counts_m,
                              colData = metadata,
                              design = ~ tissue)
dds <- DESeq(dds)
dds

res <- results(dds)
res

volcano_data <- as_tibble(res, rownames = 'gene_id')

volcano_plot <- volcano_data |> 
    ggplot(aes(x = log2FoldChange, y = -log10(padj))) +
    geom_point() +
    geom_hline(yintercept = -log10(0.05)) +
    geom_vline(xintercept = 2) +
    geom_vline(xintercept = -2) +
    theme_grey()
volcano_plot

ggsave("Notebook/plot/volcano.png", volcano_plot)

vsd <- vst(dds)
vsd

PCA <- plotPCA(vsd, intgroup = c("tissue"))

ggsave("Notebook/plot/PCA.png")