library(dplyr)
library(tidyr)
library(edgeR)

setwd("~/data/LshCas13a_RNA_cleavage/LshCas13a_M13_infection")
design_table <- read.delim("design.tsv", sep = "\t", stringsAsFactors = F, header = T)

targeting_samples <- design_table[design_table$Exp_group == "targeting", ]$Sample
nontargeting_samples <- design_table[design_table$Exp_group == "nontargeting", ]$Sample

N5E_counts_table <- read.delim("Results/Tables/Merged_ends_counts/N5E_T_vs_NT.tsv.gz", sep = "\t", header = T, stringsAsFactors = F)

N5E_counts_table$T_aveLogCPM <- aveLogCPM(N5E_counts_table[targeting_samples], prior.count = 1)
N5E_counts_table$NT_aveLogCPM <- aveLogCPM(N5E_counts_table[nontargeting_samples], prior.count = 1)
N5E_counts_table <- N5E_counts_table %>% 
  mutate(., aveLogFC = T_aveLogCPM - NT_aveLogCPM)

N5E_counts_table_CPM <- N5E_counts_table
N5E_counts_table_CPM[4:9] <- cpm(N5E_counts_table_CPM[4:9])

N5E_counts_table_CPM$T_aveCPM <- rowMeans(N5E_counts_table_CPM[targeting_samples])
N5E_counts_table_CPM$NT_aveCPM <- rowMeans(N5E_counts_table_CPM[nontargeting_samples])

std <- function(x) sd(x)/sqrt(length(x))
N5E_counts_table_CPM$TSEM <- apply(N5E_counts_table_CPM[targeting_samples], 1, std)
N5E_counts_table_CPM$NTSEM <- apply(N5E_counts_table_CPM[nontargeting_samples], 1, std)

write.table(N5E_counts_table_CPM, "Results/Tables/Merged_ends_counts/N5E_T_vs_NT_CPM_normalized.tsv", sep = "\t", quote = F, row.names = F)
