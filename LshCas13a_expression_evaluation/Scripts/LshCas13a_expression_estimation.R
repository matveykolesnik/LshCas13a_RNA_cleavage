library(edgeR)
library(dplyr)
library(tibble)
library(data.table)
library(matrixStats)

WD <- "~/data/LshCas13a_RNA_cleavage/LshCas13a_expression_evaluation/"
setwd(WD)

#load counts table
E_coli_counts_table <- "Results/Tables/HTSeqCounts/E_coli_TAP_untreated.tsv"
E_coli_counts_table.df <- read.csv(E_coli_counts_table, sep = "\t", row.names="GeneID", comment.char = "#")

sample_ids <- E_coli_counts_table.df %>% 
  select(-Name) %>% 
  colnames(.,)

gene_names <- E_coli_counts_table.df %>% 
  rownames_to_column(var = "GeneID") %>% 
  select(c("GeneID", "Name"))

E_coli_counts_table.cpm <- cpm(E_coli_counts_table.df[sample_ids]) %>% 
  as.data.frame() %>% 
  filter(rowSums(. >= 1) >= 3)

E_coli_counts_table.cpm$sd <- apply(E_coli_counts_table.cpm[sample_ids], 1, sd)
E_coli_counts_table.cpm$mean <- apply(E_coli_counts_table.cpm[sample_ids], 1, mean)
E_coli_counts_table.cpm <- E_coli_counts_table.cpm %>% 
  mutate(cv = sd/mean) %>% 
  rownames_to_column(var = "GeneID")

E_coli_counts_table.cpm <- inner_join(gene_names, E_coli_counts_table.cpm, by="GeneID") %>% 
  arrange(., cv)

low_cv_genes <- E_coli_counts_table.cpm[0:100,]
write.table(low_cv_genes, file = "Results/Tables/E_coli_low_cv_genes.tsv", sep = "\t", row.names = F, quote = F)

library(ggplot2)
library(ggrepel)

ggplot(data = E_coli_counts_table.cpm, aes(x = log10(mean), y = cv)) +
  geom_point()

reference_genes_list <- c("mrdB", "murA", "dnaE", "dnaA", "ffh")

reference_genes_data <- E_coli_counts_table.cpm %>% 
  filter(., Name %in% reference_genes_list)

E_coli_reference_genes_GM <- exp(mean(log(reference_genes_data$mean)))

E_coli_counts_table.cpm <- E_coli_counts_table.cpm %>% 
  mutate(normalizedmean = mean/E_coli_reference_genes_GM)
