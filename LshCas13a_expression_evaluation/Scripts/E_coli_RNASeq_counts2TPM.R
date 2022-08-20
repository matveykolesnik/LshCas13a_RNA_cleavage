WD <- "~/data/LshCas13a_RNA_cleavage/LshCas13a_expression_evaluation/"
setwd(WD)

library(dplyr)
library(tibble)
library(rtracklayer)
library(data.table)

#load raw fragments counts
E_coli_counts_table <- "Results/Tables/HTSeqCounts/E_coli_TAP_untreated.tsv"
E_coli_counts_table.df <- read.csv(E_coli_counts_table, sep = "\t", comment.char = "#")

#load annotation to extract genes lengths
E_coli_annotation <- "Annotations/E_coli_merged_annotation.gff3"
E_coli_annotation.df <- as.data.frame(readGFF(E_coli_annotation)) %>% 
  filter(., type=="gene") %>% 
  mutate(., length=end-start+1) %>% 
  select(GeneID = ID, Name, length)

E_coli_counts_table_length.df <- inner_join(E_coli_counts_table.df, E_coli_annotation.df, by=c("GeneID", "Name"))

#calculate TPM values
E_coli_counts_matrix <- E_coli_counts_table_length.df %>% 
  select(-c("GeneID", "Name", "length")) %>% 
  as.matrix()
rownames(E_coli_counts_matrix) <- E_coli_counts_table_length.df$GeneID

norm_by_length <- E_coli_counts_matrix/E_coli_counts_table_length.df$length
E_coli_tpm_dt <- t(t(norm_by_length)*10**6/colSums(norm_by_length)) %>% 
  as.data.frame() %>% 
  rownames_to_column(., var = "GeneID")
E_coli_tpm_dt$Name <- E_coli_counts_table_length.df$Name

#save to file
fwrite(E_coli_tpm_dt, file = "Results/Tables/E_coli_RNASeq_TPM_normalized.tsv", sep = "\t", quote = F, row.names = F)

#calculate sd and cv
sample_ids <- E_coli_tpm_dt %>% 
  select(-c("GeneID", "Name")) %>% 
  colnames(.,)

nt_samples <- c("KSIJ04", "KSIJ05", "KSIJ06")

E_coli_tpm_dt$sd <- apply(E_coli_tpm_dt[sample_ids], 1, sd)
E_coli_tpm_dt$mean <- apply(E_coli_tpm_dt[sample_ids], 1, mean)
E_coli_tpm_dt$nt_mean <- apply(E_coli_tpm_dt[nt_samples], 1, mean)
E_coli_tpm_dt <- E_coli_tpm_dt %>% 
  mutate(cv = sd/mean) %>% 
  arrange(., cv)

#select top-100 genes with low cv and export the names to file
write.table(file = "Results/Tables/E_coli_low_cv_genes_names.txt", E_coli_tpm_dt[0:100,]$Name, quote = F, row.names = F, col.names = F, sep = "")

#load the results of tblastn search for the homologs of low-cv genes in L. seeligeri genome
L_seeligeri_tblastn_search <- "Results/Tables/E_coli_low_cv_genes_tblastn_vs_L_seeligeri.tsv"
L_seeligeri_tblastn_search.df <- read.delim(L_seeligeri_tblastn_search, sep = "\t", col.names = c("qseqid",
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

#load l. seeligeri counts
L_seeligeri_counts_table <- "Results/Tables/HTSeqCounts/L_seeligeri_SLCC3954_WT_dCRISPR.tsv"
L_seeligeri_counts_table.df <- read.csv(L_seeligeri_counts_table, sep = "\t", comment.char = "#")

#load L. seeligeri annotation to get genes lengths
L_seeligeri_annotation <- "Annotations/L_seeligeri_SLCC3954.gff3"
L_seeligeri_annotation.df <- as.data.frame(readGFF(L_seeligeri_annotation)) %>% 
  filter(., type=="gene") %>% 
  mutate(., length=end-start+1) %>% 
  select(GeneID = ID, Name, length)

L_seeligeri_counts_table_length.df <- inner_join(L_seeligeri_counts_table.df, L_seeligeri_annotation.df, by=c("GeneID", "Name"))

#calculate TPM values
L_seeligeri_counts_matrix <- L_seeligeri_counts_table_length.df %>% 
  select(-c("GeneID", "Name", "length")) %>% 
  as.matrix()
rownames(L_seeligeri_counts_matrix) <- L_seeligeri_counts_table_length.df$GeneID

norm_by_length <- L_seeligeri_counts_matrix/L_seeligeri_counts_table_length.df$length
L_seeligeri_tpm_df <- t(t(L_seeligeri_counts_matrix)*10**6/colSums(L_seeligeri_counts_matrix)) %>% 
  as.data.frame() %>% 
  rownames_to_column(., var = "GeneID")
L_seeligeri_tpm_df$Name <- L_seeligeri_counts_table_length.df$Name

#attempt to check the correlation
L_seeligeri_WT <- L_seeligeri_tpm_df %>% 
  select(Name, LseWT = WT) %>% 
  mutate(logLseWT = log10(LseWT))
E_coli_mean <- E_coli_tpm_dt %>% 
  select(Name, EcoNTMean = nt_mean) %>% 
  mutate(logEcoNTMean = log10(EcoNTMean))
E_coli_mean["Name"][E_coli_mean["Name"] == "LshCas13a"] <- "cas13a"
Eco_Lse_tpm <- inner_join(E_coli_mean, L_seeligeri_WT, by="Name")
print(nrow(Eco_Lse_tpm))

library(ggplot2)
library(ggrepel)
library(ggpubr)
library(data.table)

reference_genes <- E_coli_tpm_dt[0:100,]$Name
Eco_Lse_tpm_labeled <- Eco_Lse_tpm %>% 
  mutate(ref = ifelse(Name=="cas13a", Name, ""))
  #mutate(ref = ifelse(Name %in% c("cas13a", reference_genes), Name, "")) %>% 
  #mutate(fill = ifelse(Name %in% reference_genes, "ref", "other"))


correlation_plot <- ggplot(Eco_Lse_tpm_labeled, aes(x = logEcoNTMean, y = logLseWT)) +
  #geom_point(aes(color=fill)) +
  geom_point()+
  geom_smooth(method=lm) +
  #scale_color_manual(values = c("black", "red")) +
  geom_text_repel(aes(label = ref)) +
  xlab("lgTPM of E. coli gene transcripts") +
  ylab("lgTPM of L. seeligeri gene transcripts") +
  stat_cor(method = "spearman") +
  theme_bw()

ggsave(filename = "Results/Pictures/Eco_Lse_cor_plot_spearman_test.png", plot = correlation_plot, dpi = "retina", height = 10, width = 10)
