library(DESeq2)
library(dplyr)

WD <- "~/data/LshCas13a_RNA_cleavage/LshCas13a_expression_evaluation/"
setwd(WD)

E_coli_counts_table <- "Results/Tables/HTSeqCounts/E_coli_TAP_untreated.tsv"
E_coli_counts_table_coldata <- "Results/Tables/HTSeqCounts/E_coli_TAP_untreated_colnames.tsv"

E_coli_counts_matrix <- read.csv(E_coli_counts_table, sep = "\t", row.names="Name", comment.char = "#") %>% 
  select(-GeneID) %>% 
  as.matrix()
E_coli_counts_table_coldata.df <- read.csv(E_coli_counts_table_coldata, sep = "\t", row.names = "SampleID") %>% 
  mutate(Condition = factor(Condition))

all(rownames(E_coli_counts_table_coldata.df) == colnames(E_coli_counts_matrix))

dds <- DESeqDataSetFromMatrix(countData = E_coli_counts_matrix,
                              colData = E_coli_counts_table_coldata.df,
                              design = ~ Condition)
dds <- DESeq(dds)

res <- results(dds, contrast = c("Condition", "targeting", "nontargeting"))

resOrdered <- res[order(res$pvalue),]

resOrdered

write.table(as.data.frame(resOrdered), 
            file="Results/Tables/DESeq2_results/nontargeting_vs_targeting_dge.tsv", 
            sep = "\t",
            quote = F)

res$lfcSE
