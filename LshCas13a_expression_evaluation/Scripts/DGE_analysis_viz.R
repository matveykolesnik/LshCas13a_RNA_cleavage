library(dplyr)
library(ggplot2)
library(ggrepel)
library(data.table)

WD <- "~/data/LshCas13a_RNA_cleavage/LshCas13a_expression_evaluation/"
setwd(WD)

EcC3000_DESeq_file <- "Results/Tables/DESeqRes/EcC3000_DESeqRes.tsv"
Ecd10LVM_DESeq_file <- "Results/Tables/DESeqRes/Ecd10LVM_DESeqRes.tsv"

EcC3000_DESeq.dt <- fread(EcC3000_DESeq_file)
Ecd10LVM_DESeq.dt <- fread(Ecd10LVM_DESeq_file)

E_coli_ribosomal_genes_file <- "Results/Tables/E_coli_ribosomal_proteins_.txt"
E_coli_ribosomal_genes_list <- scan(E_coli_ribosomal_genes_file, what = "character")

EcC3000_DESeq.dt.labeled <- EcC3000_DESeq.dt %>% 
  mutate(label = ifelse(Name == "rpoS", Name, "")) %>% 
  mutate(Gene = ifelse(Name %in% E_coli_ribosomal_genes_list, "Ribosomal\nprotein genes",
                       ifelse(Name == "rpoS", "rpoS", "Other"))) %>% 
  arrange(Gene)

color_scheme <- c("Ribosomal\nprotein genes" = "red", "Other" = "grey", "rpoS" = "black")

EcC3000_volcano <- ggplot(data = EcC3000_DESeq.dt.labeled, aes(x = log2FoldChange, y = -log10(padj), color=Gene)) +
  geom_point() +
  scale_color_manual(values = color_scheme) +
  geom_text_repel(aes(label = label), max.overlaps = 10000, show.legend = F) +
  xlab("log2FC") +
  ylab("-log10(FDR)") +
  theme_bw()

ggsave("Results/Pictures/EcC3000_DESeq2_volcano_rp_genes_rpoS_colored.png", plot = EcC3000_volcano, height = 5, width = 5)

Type_I_antitoxin_rnas <- c("rdlA", "rdlB", "rdlC", "rdlD", "istR", "symR", "sibA", "sibB", "sibC", "sibD", "sibE", "ohsC", "ralA", "agrA")
color_scheme1 <- c("antitoxin RNA" = "red", "Other" = "grey")

EcC3000_DESeq.dt.labeled_TAS <- EcC3000_DESeq.dt %>% 
  mutate(label = ifelse(Name %in% Type_I_antitoxin_rnas, "antitoxin RNA", "Other")) %>% 
  arrange(desc(label))

EcC3000_volcano_Type_I_TAS <- ggplot(data = EcC3000_DESeq.dt.labeled_TAS, aes(x = log2FoldChange, y = -log10(padj), color=label)) +
  geom_point() +
  scale_color_manual(values = color_scheme1) +
  xlab("log2FC") +
  ylab("-log10(FDR)") +
  theme_bw()

ggsave("Results/Pictures/EcC3000_DESeq2_volcano_Type_I_TAS_colored.png", plot = EcC3000_volcano_Type_I_TAS, height = 5, width = 5)

Ecd10LVM_DESeq.dt.labeled <- Ecd10LVM_DESeq.dt %>% 
  mutate(label = ifelse(Name == "rtcB", Name, "")) %>% 
  mutate(Gene = ifelse(Name == "rtcB", Name, "Other")) %>% 
  arrange(label)

color_scheme2 <- c("rtcB" = "red", "Other" = "grey")

Ecd10LVM_volcano <- ggplot(data = Ecd10LVM_DESeq.dt.labeled, aes(x = log2FoldChange, y = -log10(padj), color=Gene)) +
  geom_point() +
  scale_color_manual(values = color_scheme2) +
  xlab("log2FC") +
  ylab("-log10(FDR)") +
  theme_bw()

ggsave("Results/Pictures/Ecd10LVM_DESeq2_volcano_rtcB_colored.png", plot = Ecd10LVM_volcano, height = 5, width = 5)

Gourse_2019_down_genes <- c("ydgI")
Gourse_2019_up_genes <- c("argC", "argB", "argH", "argA", "argE", "argF", "argI")

EcC3000_DESeq.dt.Gourse <- EcC3000_DESeq.dt %>% 
  mutate(Gene = ifelse(Name %in% Gourse_2019_down_genes, "downregulated (Sanchez-Vazquez et al., 2019)",
                       ifelse(Name %in% Gourse_2019_up_genes, "upregulated (Sanchez-Vazquez et al., 2019)", "Other"))) %>% 
  arrange(Gene)
color_scheme3 <- c("downregulated (Sanchez-Vazquez et al., 2019)" = "red",
                   "upregulated (Sanchez-Vazquez et al., 2019)" = "orange",
                   "Other" = "grey")


EcC3000_volcano_Goerse_data <- ggplot(data = EcC3000_DESeq.dt.Gourse, aes(x = log2FoldChange, y = -log10(padj), color=Gene)) +
  geom_point() +
  scale_color_manual(values = color_scheme3) +
  xlab("log2FC") +
  ylab("-log10(FDR)") +
  theme_bw()

ggsave("Results/Pictures/EcC3000_DESeq2_Goerse_data_arg.png", plot = EcC3000_volcano_Goerse_data, height = 10, width = 10)
