library(dplyr)
library(ggplot2)
library(ggrepel)

WD <- "~/data/LshCas13a_RNA_cleavage/LshCas13a_expression_evaluation/"
setwd(WD)

N=100

DESeq_resuls <- "Results/Tables/DESeq2_results/nontargeting_vs_targeting_dge.tsv"
DESeq_resuls.df <- read.delim(DESeq_resuls, sep = "\t", header = T, stringsAsFactors = F) %>% 
  mutate(Name = c(rownames(.,)[1:N], rep("", nrow(.)-N)))

ggplot(DESeq_resuls.df, aes(x = log2FoldChange, y = -log10(padj))) +
  geom_point() +
  geom_text_repel(aes(label=Name))

DESeq_resuls_poslfc.df <- DESeq_resuls.df %>% 
  filter(log2FoldChange > 0)
