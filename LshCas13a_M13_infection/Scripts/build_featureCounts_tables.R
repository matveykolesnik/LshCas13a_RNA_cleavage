setwd("~/data/LshCas13a_RNA_cleavage/LshCas13a_M13_infection")
library(Rsubread)

EcC3000M13BAMFilesDir <- "Alignments/BAM"
EcC3000M13BAMFilesList <- sapply(list.files(EcC3000M13BAMFilesDir, pattern = "\\.bam$"), 
                                 function(x) paste0(EcC3000M13BAMFilesDir, "/", x))

AnnotationFile <- "Annotations/Summarized_annotation.gff"

EcC3000M13_featureCounts <- featureCounts(files = EcC3000M13BAMFilesList, 
                                          annot.ext = AnnotationFile, 
                                          isGTFAnnotationFile = T,
                                          GTF.featureType = "gene",
                                          GTF.attrType = "ID",
                                          GTF.attrType.extra = "Name",
                                          allowMultiOverlap = T,
                                          strandSpecific = 0,
                                          isPairedEnd = F,
                                          # reportReads = "CORE",
                                          # reportReadsPath = "Results/featureCounts_CORE_files",
                                          nthreads = 16)

library(data.table)
library(dplyr)

EcC3000M13_featureCounts_raw <- left_join(as.data.table(EcC3000M13_featureCounts$counts, keep.rownames = "GeneID"),
                                          EcC3000M13_featureCounts$annotation,
                                          by="GeneID") %>% 
  select(-c("Start", "End", "Strand", "Length", "Chr")) %>% 
  select(GeneID, Name, everything())
names(EcC3000M13_featureCounts_raw) <- sub("_S[0-9]+_sorted.bam", "", names(EcC3000M13_featureCounts_raw))
fwrite(x = EcC3000M13_featureCounts_raw, file = "Results/Tables/featureCountsTables/EcC3000M13_featureCounts_raw_counts.tsv",
       sep = "\t", quote = F, row.names = F)

#build TPM tables
EcC3000M13_featureCounts_matrix <- as.matrix(EcC3000M13_featureCounts$counts)
EcC3000M13_featureCounts_matrix.norm_by_length <- EcC3000M13_featureCounts_matrix/EcC3000M13_featureCounts$annotation$Length
EcC3000M13_featureCounts.TPM <- t(t(EcC3000M13_featureCounts_matrix.norm_by_length)*10**6/colSums(EcC3000M13_featureCounts_matrix.norm_by_length)) %>% 
  as.data.table(., keep.rownames = "GeneID") %>% 
  left_join(., EcC3000M13_featureCounts$annotation, by="GeneID") %>% 
  select(-c("Start", "End", "Strand", "Length", "Chr")) %>% 
  select(GeneID, Name, everything())
names(EcC3000M13_featureCounts.TPM) <- sub("_S[[:digit:]]+_sorted.bam", "", names(EcC3000M13_featureCounts.TPM))
fwrite(x = EcC3000M13_featureCounts.TPM, file = "Results/Tables/featureCountsTables/EcC3000M13_featureCounts_TPM_normalized.tsv",
       sep = "\t", quote = F, row.names = F)

T_samples_list <- c("M13_T1", "M13_T2", "M13_T3")
NT_samples_list <- c("M13_NT1", "M13_NT2", "M13_NT3")

EcC3000M13_featureCounts.TPM$T_means <- apply(EcC3000M13_featureCounts.TPM[, ..T_samples_list], 1, mean)
EcC3000M13_featureCounts.TPM$NT_means <- apply(EcC3000M13_featureCounts.TPM[, ..NT_samples_list], 1, mean)

M13_genes <- EcC3000M13_featureCounts.TPM[startsWith(EcC3000M13_featureCounts.TPM$Name, "M13")]$Name

EcC3000M13_featureCounts.TPM_means <- EcC3000M13_featureCounts.TPM %>% 
  select(., c(Name, T_means, NT_means)) %>% 
  filter_at(., vars(T_means, NT_means), all_vars((.) != 0)) %>% 
  mutate(Gene = ifelse(Name %in% M13_genes, "M13 genes", "Other")) %>% 
  arrange(desc(Gene))

library(ggplot2)
library(ggrepel)
library(ggpubr)

color_scheme <- c("M13 genes" = "red", "Other" = "grey")

EcC3000_M13_T_vs_NT <- ggplot(EcC3000M13_featureCounts.TPM_means, aes(x = log10(T_means), y = log10(NT_means), color=Gene)) +
  #geom_point(aes(color=fill)) +
  geom_point()+
  geom_abline(intercept = 0, slope = 1, color = "blue") +
  scale_color_manual(values = color_scheme) +
  xlab("lgTPM of transcripts in targeting samples") +
  ylab("lgTPM of transcripts in nontargeting samples") +
  theme_bw()

ggsave(filename = "Results/Pictures/EcC3000_M13_T_vs_NT_plot.png", plot = EcC3000_M13_T_vs_NT, width = 5, height = 5)
