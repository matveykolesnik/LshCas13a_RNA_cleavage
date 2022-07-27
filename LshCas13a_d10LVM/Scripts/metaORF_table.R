library(tidyverse)
library(ggrepel)
library(data.table)

setwd("LshCas13a_d10LVM/")

TCS_features_table <- read.delim("Results/Tables/TCS_detection_tables/LRTest_table_and_genome_features.tsv", sep = "\t", header = T, stringsAsFactors = F) %>% 
  arrange(., PValue.adj)

#select TCS located within CDS with logFC > 0, i. e. positons where 5'ends of fragments are enriched in targeting samples comparing with nontargeting samples and located in the same strand as ORF
CDS_located_TCS_poslogfc <- TCS_features_table %>% 
  filter(., MatchedFeatureType == "CDS", Strand == MatchedFeatureStrand, logFC>0) %>% 
  mutate(., DistFromORFStart = ifelse(MatchedFeatureStrand == "+", Pos-MatchedFeatureStart, MatchedFeatureEnd-Pos)) %>% #calculate distance
  filter(., row_number() <= 1000) %>% 
  rowwise() %>% 
  mutate(., DistFromORFStart_rnd = sample(0:(MatchedFeatureEnd-MatchedFeatureStart+1), 1))

limit=20

cleavage_pos <- c("-1/1")
for (i in 1:limit) {
  #print(paste0(i, "/", i+1))
  cleavage_pos <- c(cleavage_pos, paste0(i, "/", i+1))
}

barplot_df <- data.table(Dist = seq(0, limit, by=1))

barplot_df$DistFromORFReal <- apply(barplot_df, 1, function(x) nrow(CDS_located_TCS_poslogfc[CDS_located_TCS_poslogfc$DistFromORFStart == x["Dist"],]))
barplot_df$DistFromORFRnd <- apply(barplot_df, 1, function(x) nrow(CDS_located_TCS_poslogfc[CDS_located_TCS_poslogfc$DistFromORFStart_rnd == x["Dist"],]))
barplot_df$Cleavage_pos <- factor(cleavage_pos, levels = cleavage_pos)

write.table(barplot_df, "Results/Tables/Picture_source_tables/Ecd10LVM_TCS_across_metaORF.tsv", row.names = F, quote = F, sep = "\t")

top_1000_limit = as.integer(max(apply(CDS_located_TCS_poslogfc, 1, function(x) max(c(x["DistFromORFStart"], x["DistFromORFStart_rnd"])))))

barplot_df_top1000 <- data.table(Dist = seq(0, top_1000_limit, by=1))
barplot_df_top1000$DistFromORFReal <- apply(barplot_df_top1000, 1, function(x) nrow(CDS_located_TCS_poslogfc[CDS_located_TCS_poslogfc$DistFromORFStart == x["Dist"],]))
barplot_df_top1000$DistFromORFRnd <- apply(barplot_df_top1000, 1, function(x) nrow(CDS_located_TCS_poslogfc[CDS_located_TCS_poslogfc$DistFromORFStart_rnd == x["Dist"],]))

write.table(barplot_df_top1000, "Results/Tables/Picture_source_tables/Ecd10LVM_TCS_across_metaORF_big_scale.tsv", row.names = F, quote = F, sep = "\t")
