setwd("~/data/LshCas13a_RNA_cleavage/LshCas13a_M13_infection")

library(dplyr)
library(ggplot2)
library(ggrepel)
library(data.table)

EcC3000M13_DESeq2_file <- "Results/Tables/DESeqRes/EcC3000M13_DESeqRes.tsv"
EcC3000M13_DESeq2.dt <- fread(EcC3000M13_DESeq2_file)

M13_genes <- EcC3000M13_DESeq2.dt[startsWith(EcC3000M13_DESeq2.dt$Name, "M13")]$Name

EcC3000M13_DESeq2.dt.labeled <- EcC3000M13_DESeq2.dt %>% 
  mutate(Gene = ifelse(Name %in% M13_genes, "M13 genes", "Other")) %>% 
  arrange(desc(Gene))

color_scheme <- c("M13 genes" = "red", "Other" = "grey")


EcC3000M13_DESeq2_volcano <- ggplot(data = EcC3000M13_DESeq2.dt.labeled, aes(x = log2FoldChange, y = -log10(padj), color=Gene)) +
  geom_point() +
  scale_color_manual(values = color_scheme) +
  xlab("log2FC") +
  ylab("-log10(FDR)") +
  theme_bw()

ggsave("Results/Pictures/EcC3000M13_M13_genes.png", plot = EcC3000M13_DESeq2_volcano, height = 5, width = 5)
