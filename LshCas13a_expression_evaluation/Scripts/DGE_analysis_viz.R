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
Type_II_antitoxin_genes <- fread("MG1655_type_II_antitoxins.txt")$antitoxin

color_scheme1 <- c("Type I antitoxin RNA" = "red", 
                   "Type II antitoxin mRNA" = "yellow",
                   "Other" = "grey")

EcC3000_DESeq.dt.labeled_TAS <- EcC3000_DESeq.dt %>% 
  mutate(label = ifelse(Name %in% Type_I_antitoxin_rnas, "Type I antitoxin RNA", 
                        ifelse(Name %in% Type_II_antitoxin_genes, "Type II antitoxin mRNA", "Other"))) %>% 
  arrange(label)

EcC3000_volcano_Type_I_TAS <- ggplot(data = EcC3000_DESeq.dt.labeled_TAS, aes(x = log2FoldChange, y = -log10(padj), color=label)) +
  geom_point() +
  scale_color_manual(values = color_scheme1) +
  xlab("log2FC") +
  ylab("-log10(FDR)") +
  geom_hline(yintercept=-log10(0.001), linetype="dashed", color = "black") +
  theme_bw() +
  theme(legend.position="none")

ggsave(filename = "Results/Pictures/Type_I_antitoxin_RNAs_expression/EcC3000_DESeq2_volcano_Type_I_and_II_TAS_colored.png", 
       plot = EcC3000_volcano_Type_I_TAS, height = 8, width = 8)

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

#build separate plots
EcC3000_DESeq.dt.Gourse_upregulated <- EcC3000_DESeq.dt %>% 
  mutate(Gene = ifelse(Name %in% Gourse_2019_up_genes, "upregulated \n(Sanchez-Vazquez et al., 2019)", "Other")) %>% 
  arrange(Gene)
color_scheme_up = c("upregulated \n(Sanchez-Vazquez et al., 2019)" = "red", "other" = "grey")

EcC3000_DESeq.dt.Gourse_upregulated_plot <- ggplot(data = EcC3000_DESeq.dt.Gourse_upregulated, aes(x = log2FoldChange, y = -log10(padj), color=Gene)) +
  geom_point() +
  scale_color_manual(values = color_scheme_up) +
  xlab("log2FC") +
  ylab("-log10(FDR)") +
  geom_hline(yintercept = -log10(0.001), linetype="dashed", colour="black") +
  theme_bw() +
  theme(text = element_text(size = 20))

ggsave("Results/Pictures/Stress_response/Goerse_et_al_upregulated_genes.png", 
       plot = EcC3000_DESeq.dt.Gourse_upregulated_plot, height = 8, width = 8)

#downregulated
EcC3000_DESeq.dt.Gourse_downregulated <- EcC3000_DESeq.dt %>% 
  mutate(Gene = ifelse(Name %in% Gourse_2019_down_genes, "downregulated \n(Sanchez-Vazquez et al., 2019)", "Other")) %>% 
  arrange(desc(Gene))
color_scheme_down = c("downregulated \n(Sanchez-Vazquez et al., 2019)" = "red", "other" = "grey")

EcC3000_DESeq.dt.Gourse_downregulated_plot <- ggplot(data = EcC3000_DESeq.dt.Gourse_downregulated, aes(x = log2FoldChange, y = -log10(padj), color=Gene)) +
  geom_point() +
  scale_color_manual(values = color_scheme_down) +
  xlab("log2FC") +
  ylab("-log10(FDR)") +
  geom_hline(yintercept = -log10(0.001), linetype="dashed", colour="black") +
  theme_bw() +
  theme(text = element_text(size = 20))

ggsave("Results/Pictures/Stress_response/Goerse_et_al_downregulated_genes.png",
       plot = EcC3000_DESeq.dt.Gourse_downregulated_plot, height = 8, width = 8)

library(rtracklayer)

Annotation <- "Annotations/E_coli_merged_annotation.gff3"

Annotation.df <- as.data.frame(readGFF(Annotation)) %>% 
  filter(type == "gene") %>% 
  mutate(gene_biotype = sub("protein_coding", "CDS", gene_biotype)) %>%
  mutate(gene_biotype = ifelse(gene_biotype %in% c("ncRNA", NA, "rRNA"), "Other", gene_biotype)) %>% 
  select(c("Name", "gene_biotype")) %>% 
  mutate(gene_biotype = factor(gene_biotype, levels = c("tRNA", "CDS", "Other")))

Ecd10LVM_DESeq.dt.gene_biotypes <- inner_join(Ecd10LVM_DESeq.dt, Annotation.df)

fwrite(Ecd10LVM_DESeq.dt.gene_biotypes, file = "Results/Tables/Ecd10LVM_DESeq.dt.gene_biotypes.tsv", sep = "\t", row.names = F)

Ecd10LVM_DESeq_volcano <- ggplot(data = Ecd10LVM_DESeq.dt.gene_biotypes, aes(x = log2FoldChange, y = -log10(padj), color=gene_biotype)) +
  geom_point() +
  scale_color_manual(values = c("tRNA" = "red", "CDS" = "seagreen", "Other" = "grey"), name = "") +
  xlab("log2FC") +
  ylab("-log10(FDR)") +
  geom_hline(yintercept = -log10(0.001), linetype="dashed", colour="black") +
  theme_bw() +
  theme(text = element_text(size = 20))

ggsave("Results/Pictures/Ecd10LVM_DESeq2_volcano_tRNA_CDS_colored.png", Ecd10LVM_DESeq_volcano, height = 8, width = 8)

