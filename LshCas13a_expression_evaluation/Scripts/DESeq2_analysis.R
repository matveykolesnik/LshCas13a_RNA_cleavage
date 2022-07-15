setwd("LshCas13a_expression_evaluation/")
library(DESeq2)
library(dplyr)
library(data.table)

EcC3000_raw_counts_file <- "Results/Tables/featureCountsTables/EcC3000_featureCounts_raw_counts.tsv"
Ecd10LVM_raw_counts_file <- "Results/Tables/featureCountsTables/Ecd10LVM_featureCounts_raw_counts.tsv"

EcC3000_conditions_file <- "Results/Tables/EcC3000_conditions.tsv"
Ecd10LVM_conditions_file <- "Results/Tables/Ecd10LVM_conditions.tsv"

RunDESeq <- function(CountsFile, ConditionsFile) {
  counts_matrix <- read.csv(CountsFile, sep = "\t", row.names="Name", comment.char = "#") %>% 
    select(-GeneID) %>% 
    as.matrix()
  conditions <- read.csv(ConditionsFile, sep = "\t", row.names = "SampleID") %>% 
    mutate(Condition = factor(Condition))
  if (!all(rownames(conditions) == colnames(counts_matrix))) {
    print("Rownames of the conditions table do not correspond to colnames of the count table")
    return(F)
  }
  dds <- DESeqDataSetFromMatrix(countData = counts_matrix,
                                colData = conditions,
                                design = ~ Condition)
  dds <- DESeq(dds)
  res <- results(dds, contrast = c("Condition", "targeting", "nontargeting"))
  resOrdered <- res[order(res$pvalue),]
  return(as.data.frame(resOrdered))
}

EcC3000_DESeqRes <- RunDESeq(CountsFile = EcC3000_raw_counts_file, ConditionsFile = EcC3000_conditions_file) %>% 
  as.data.table(., keep.rownames = "Name")
Ecd10LVM_DESeqRes <- RunDESeq(CountsFile = Ecd10LVM_raw_counts_file, ConditionsFile = Ecd10LVM_conditions_file) %>% 
  as.data.table(., keep.rownames = "Name")

fwrite(x = EcC3000_DESeqRes, file = "Results/Tables/DESeqRes/EcC3000_DESeqRes.tsv", sep = "\t", quote = F, row.names = F)
fwrite(x = Ecd10LVM_DESeqRes, file = "Results/Tables/DESeqRes/Ecd10LVM_DESeqRes.tsv", sep = "\t", quote = F, row.names = F)
