WD <- "~/data/LshCas13a_RNA_cleavage/LshCas13a_expression_evaluation/"
setwd(WD)

library(dplyr)

E_coli_low_cv_genes_homologs <- "Results/Tables/E_coli_low_cv_genes_tblastn_vs_L_shahii.tsv"
E_coli_low_cv_genes <- "Results/Tables/E_coli_low_cv_genes.tsv"

E_coli_low_cv_genes.df <- read.delim(E_coli_low_cv_genes) %>% 
  select(c("Name", "cv"))

E_coli_low_cv_genes_homologs.df <- read.delim(E_coli_low_cv_genes_homologs, 
                                              col.names = c("qseqid",
                                                            "sseqid",
                                                            "pident",
                                                            "length",
                                                            "mismatch",
                                                            "gapopen",
                                                            "qstart",
                                                            "qend",
                                                            "sstart",
                                                            "send",
                                                            "evalue",
                                                            "bitscore")) %>% 
  arrange(., evalue)

E_coli_low_cv_genes_homologs_cv.df <- inner_join(E_coli_low_cv_genes_homologs.df, E_coli_low_cv_genes.df, by=c("qseqid" = "Name")) %>% 
  filter(evalue <= 1e-9)
