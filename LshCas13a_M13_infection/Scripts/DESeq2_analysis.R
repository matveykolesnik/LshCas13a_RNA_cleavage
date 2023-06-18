setwd("~/data/LshCas13a_RNA_cleavage/LshCas13a_M13_infection/")
library(DESeq2)
library(dplyr)
library(data.table)

EcC3000M13_raw_counts_file <- "Results/Tables/featureCountsTables/EcC3000M13_featureCounts_raw_counts.tsv"
EcC3000M13_conditions_file <- "Results/Tables/featureCountsTables/EcC3000M13_conditions.txt"

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

EcC3000M13M13_DESeqRes <- RunDESeq(CountsFile = EcC3000M13_raw_counts_file, ConditionsFile = EcC3000M13_conditions_file) %>% 
  as.data.table(., keep.rownames = "Name")

fwrite(x = EcC3000M13M13_DESeqRes, file = "Results/Tables/DESeqRes/EcC3000M13_DESeqRes.tsv", sep = "\t", quote = F, row.names = F)
