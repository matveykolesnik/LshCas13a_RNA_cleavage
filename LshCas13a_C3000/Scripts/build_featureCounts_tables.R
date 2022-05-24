setwd("LshCas13a_expression_evaluation")
library(Rsubread)

EcC3000BAMFilesDir <- "Alignments/E_coli_C3000_BAM"
Ecd10LVMBAMFilesDir <- "Alignments/E_coli_d10LVM_BAM"

EcC0000BAMFilesList <- sapply(list.files(EcC3000BAMFilesDir, pattern = "\\.bam$"), 
                              function(x) paste0(EcC3000BAMFilesDir, "/", x))
Ecd10LVMBAMFilesList <- sapply(list.files(Ecd10LVMBAMFilesDir, pattern = "\\.bam$"), 
                               function(x) paste0(Ecd10LVMBAMFilesDir, "/", x))

AnnotationFile <- "Annotations/E_coli_merged_annotation.gff3"

EcC3000_featureCounts <- featureCounts(files = EcC0000BAMFilesList, 
                                       annot.ext = AnnotationFile, 
                                       isGTFAnnotationFile = T,
                                       GTF.featureType = "gene",
                                       GTF.attrType = "ID",
                                       GTF.attrType.extra = "Name",
                                       allowMultiOverlap = T,
                                       strandSpecific = 1,
                                       isPairedEnd = T,
                                       countReadPairs = T,
                                       requireBothEndsMapped = T,
                                       reportReads = "CORE",
                                       reportReadsPath = "Results/featureCounts_CORE_files",
                                       nthreads = 16)

Ecd10LVM_featureCounts <- featureCounts(files = Ecd10LVMBAMFilesList, 
                                        annot.ext = AnnotationFile, 
                                        isGTFAnnotationFile = T,
                                        GTF.featureType = "gene",
                                        GTF.attrType = "ID",
                                        GTF.attrType.extra = "Name",
                                        allowMultiOverlap = T,
                                        strandSpecific = 1,
                                        isPairedEnd = T,
                                        countReadPairs = T,
                                        requireBothEndsMapped = T,
                                        reportReads = "CORE",
                                        reportReadsPath = "Results/featureCounts_CORE_files",
                                        nthreads = 16)

library(edgeR)
library(data.table)
library(dplyr)

EcC3000_featureCounts_raw <- left_join(as.data.table(EcC3000_featureCounts$counts, keep.rownames = "GeneID"),
                                       EcC3000_featureCounts$annotation,
                                       by="GeneID") %>% 
  select(-c("Start", "End", "Strand", "Length", "Chr")) %>% 
  select(GeneID, Name, everything())
names(EcC3000_featureCounts_raw) <- sub("_sorted.bam", "", names(EcC3000_featureCounts_raw))
fwrite(x = EcC3000_featureCounts_raw, file = "Results/Tables/featureCountsTables/EcC3000_featureCounts_raw_counts.tsv",
       sep = "\t", quote = F, row.names = F)

EcC3000_featureCounts_matrix <- as.matrix(EcC3000_featureCounts$counts)
EcC3000_featureCounts_matrix.norm_by_length <- EcC3000_featureCounts_matrix/EcC3000_featureCounts$annotation$Length
EcC3000_featureCounts.TPM <- t(t(EcC3000_featureCounts_matrix.norm_by_length)*10**6/colSums(EcC3000_featureCounts_matrix.norm_by_length)) %>% 
  as.data.table(., keep.rownames = "GeneID") %>% 
  left_join(., EcC3000_featureCounts$annotation, by="GeneID") %>% 
  select(-c("Start", "End", "Strand", "Length", "Chr")) %>% 
  select(GeneID, Name, everything())
names(EcC3000_featureCounts.TPM) <- sub("_sorted.bam", "", names(EcC3000_featureCounts.TPM))
fwrite(x = EcC3000_featureCounts.TPM, file = "Results/Tables/featureCountsTables/EcC3000_featureCounts_TPM_normalized.tsv",
       sep = "\t", quote = F, row.names = F)

#the same procedure with d10LVM samples
Ecd10LVM_featureCounts_raw <- left_join(as.data.table(Ecd10LVM_featureCounts$counts, keep.rownames = "GeneID"),
                                       Ecd10LVM_featureCounts$annotation,
                                       by="GeneID") %>% 
  select(-c("Start", "End", "Strand", "Length", "Chr")) %>% 
  select(GeneID, Name, everything())
names(Ecd10LVM_featureCounts_raw) <- sub("_sorted.bam", "", names(Ecd10LVM_featureCounts_raw))
fwrite(x = Ecd10LVM_featureCounts_raw, file = "Results/Tables/featureCountsTables/Ecd10LVM_featureCounts_raw_counts.tsv",
       sep = "\t", quote = F, row.names = F)

Ecd10LVM_featureCounts_matrix <- as.matrix(Ecd10LVM_featureCounts$counts)
Ecd10LVM_featureCounts_matrix.norm_by_length <- Ecd10LVM_featureCounts_matrix/Ecd10LVM_featureCounts$annotation$Length
Ecd10LVM_featureCounts.TPM <- t(t(Ecd10LVM_featureCounts_matrix.norm_by_length)*10**6/colSums(Ecd10LVM_featureCounts_matrix.norm_by_length)) %>% 
  as.data.table(., keep.rownames = "GeneID") %>% 
  left_join(., Ecd10LVM_featureCounts$annotation, by="GeneID") %>% 
  select(-c("Start", "End", "Strand", "Length", "Chr")) %>% 
  select(GeneID, Name, everything())
names(Ecd10LVM_featureCounts.TPM) <- sub("_sorted.bam", "", names(Ecd10LVM_featureCounts.TPM))
