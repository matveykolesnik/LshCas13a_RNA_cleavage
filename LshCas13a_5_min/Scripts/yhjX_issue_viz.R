library(dplyr)
library(ggplot2)
library(ggrepel)
library(data.table)

setwd("~/data/LshCas13a_RNA_cleavage/LshCas13a_5_min")

EcC3000_5_min_DESeq2_file <- "Results/Tables/DESeqRes/EcC3000_5_min_DESeqRes.tsv"
EcC3000_5_min_DESeq2.dt <- fread(EcC3000_5_min_DESeq2_file)

downregulated_genes <- c("glpK", "glpQ", "yhjX", "treB")

EcC3000_5_min_DESeq2.dt.labeled <- EcC3000_5_min_DESeq2.dt %>% 
  mutate(label = ifelse(Name %in% downregulated_genes, Name, "")) %>% 
  arrange(label)

EcC3000_5_min_volcano_yhjX <- ggplot(data = EcC3000_5_min_DESeq2.dt.labeled, aes(x = log2FoldChange, y = -log10(padj))) +
  geom_point() +
  geom_text_repel(aes(label = label), show.legend = F) +
  xlab("log2FC") +
  ylab("-log10(FDR)") +
  theme_bw()

ggsave("Results/Pictures/EcC3000_5_min_volcano_downregulated_genes.png", plot = EcC3000_5_min_volcano_yhjX, height = 5, width = 5)
