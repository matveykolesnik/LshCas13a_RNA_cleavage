library(dplyr)
library(data.table)

setwd("~/data/LshCas13a_RNA_cleavage/tRNA_cleavage_analysis/")

#LRTable <- "../LshCas13a_C3000/Results/Tables/TCS_detection_tables/LRTest_table_and_genome_features.tsv"
LRTable <- "../LshCas13a_d10LVM/Results/Tables/TCS_detection_tables/LRTest_table_and_genome_features.tsv"

LRTable.DT <- fread(LRTable, stringsAsFactors = F) %>% 
  filter(MatchedFeatureType=="tRNA")

LRTable.DT.RelPos <- LRTable.DT %>% 
  mutate(RelPos = ifelse(Strand == "+", Pos-MatchedFeatureStart, MatchedFeatureEnd-Pos))

for (gene in unique(LRTable.DT.RelPos$MatchedFeatureGene)) {
  LRTable.DT.RelPos.subset <- LRTable.DT.RelPos %>% 
    filter(MatchedFeatureGene == gene) %>% 
    select(c(SeqID, RelPos, Strand, logFC, PValue.adj)) %>% 
    filter(logFC > 0) %>% 
    arrange(RelPos)
  fwrite(LRTable.DT.RelPos.subset, file = paste0("Results/d10LVM_LRTable_tRNA_subsets/", gene, ".tsv"), sep = "\t", row.names = F, quote = F)
}
