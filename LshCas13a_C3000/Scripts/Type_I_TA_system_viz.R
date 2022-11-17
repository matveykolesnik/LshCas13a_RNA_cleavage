library(dplyr)
library(ggplot2)
library(ggrepel)
library(data.table)

WD <- "~/data/LshCas13a_RNA_cleavage/LshCas13a_C3000"
setwd(WD)

TCS_annotated_table <- "Results/Tables/TCS_detection_tables/LRTest_table_and_genome_features.tsv"
TCS_annotated_table.dt <- fread(TCS_annotated_table)

Type_I_antitoxin_rnas <- c("rdlA", "rdlB", "rdlC", "rdlD", "istR", "symR", "sibA", "sibB", "sibC", "sibD", "sibE", "ohsC", "ralA", "agrA")

TCS_annotated_table.dt_labeled <- TCS_annotated_table.dt %>% 
  mutate(label = ifelse(MatchedFeatureGene %in% Type_I_antitoxin_rnas, "antitoxin RNA", "Other")) %>% 
  arrange(desc(label))

color_scheme <- c("antitoxin RNA" = "red", "Other" = "grey")

Type_I_antitoxin_rnas_cleavage_volcano <- ggplot(data = TCS_annotated_table.dt_labeled, aes(x = logFC, y = -log10(PValue.adj), color=label)) +
  geom_point() +
  scale_color_manual(values = color_scheme) +
  xlab("log2FC") +
  ylab("-log10(FDR)") +
  theme_bw()

ggsave("Results/Pictures/EcC3000_TCS_Type_I_antitoxin_rnas_colored.png", plot = Type_I_antitoxin_rnas_cleavage_volcano, height = 5, width = 5)
