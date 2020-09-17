library(tidyverse)
library(ggrepel)
library(data.table)

setwd("~/data/LshCas13a_RNA_cleavage/LshCas13a_in_vitro_total_RNA/")

TCS_features_table <- read.delim("Results/Tables/TCS_detection_tables/LRTest_table_and_genome_features.tsv", sep = "\t", header = T, stringsAsFactors = F) %>% 
  arrange(., PValue.adj) %>% 
  filter(., !(MatchedFeatureType %in% c("tRNA")))

TCS_features_table$feature_label <- as.character("")
t = 20
TCS_features_table[1:t,]$feature_label <- TCS_features_table[1:t,]$MatchedFeatureGene

#draw volcano plot
TCS_volcano <- ggplot(TCS_features_table, aes(x=logFC, y=-log10(PValue.adj), colour = MatchedFeatureType)) +
  geom_point() +
  geom_text_repel(aes(label=feature_label)) +
  ggtitle("d10LVM TCS volcano plot") +
  xlab(label = "log2 fold change between targeting and nontargeting samples") +
  ylab(label = "-log10 adjusted p-value") +
  guides(color=guide_legend("Type of genomic feature"))+
  theme_classic()

ggsave("hairpin_cleavage/Pictures/LshCas13a_in_vitro_total_RNA_volcano_nontRNA.png", plot = TCS_volcano, device = "png", dpi=500, scale = 1.5)
