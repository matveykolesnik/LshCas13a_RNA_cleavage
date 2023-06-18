setwd("~/data/LshCas13a_RNA_cleavage/LshCas13a_5_min")
library(Rsubread)

EcC3000_5_min_BAMFilesDir <- "Alignments/BAM"
EcC3000_5_min_BAMFilesList <- sapply(list.files(EcC3000_5_min_BAMFilesDir, pattern = "\\.bam$"), 
                                 function(x) paste0(EcC3000_5_min_BAMFilesDir, "/", x))

AnnotationFile <- "Annotations/Summarized_annotation.gff3"

EcC3000_5_min_featureCounts <- featureCounts(files = EcC3000_5_min_BAMFilesList, 
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

EcC3000_5_min_featureCounts_raw <- left_join(as.data.table(EcC3000_5_min_featureCounts$counts, keep.rownames = "GeneID"),
                                          EcC3000_5_min_featureCounts$annotation,
                                          by="GeneID") %>% 
  select(-c("Start", "End", "Strand", "Length", "Chr")) %>% 
  select(GeneID, Name, everything())

names(EcC3000_5_min_featureCounts_raw) <- sub("_S[0-9]+_sorted.bam", "", names(EcC3000_5_min_featureCounts_raw))
fwrite(x = EcC3000_5_min_featureCounts_raw, file = "Results/Tables/featureCountsTables/EcC3000_5_min_featureCounts_raw_counts.tsv",
       sep = "\t", quote = F, row.names = F)

#build TPM tables
EcC3000_5_min_featureCounts_matrix <- as.matrix(EcC3000_5_min_featureCounts$counts)
EcC3000_5_min_featureCounts_matrix.norm_by_length <- EcC3000_5_min_featureCounts_matrix/EcC3000_5_min_featureCounts$annotation$Length
EcC3000_5_min_featureCounts.TPM <- t(t(EcC3000_5_min_featureCounts_matrix.norm_by_length)*10**6/colSums(EcC3000_5_min_featureCounts_matrix.norm_by_length)) %>% 
  as.data.table(., keep.rownames = "GeneID") %>% 
  left_join(., EcC3000_5_min_featureCounts$annotation, by="GeneID") %>% 
  select(-c("Start", "End", "Strand", "Length", "Chr")) %>% 
  select(GeneID, Name, everything())

names(EcC3000_5_min_featureCounts.TPM) <- sub("_S[[:digit:]]+_sorted.bam", "", names(EcC3000_5_min_featureCounts.TPM))
fwrite(x = EcC3000_5_min_featureCounts.TPM, file = "Results/Tables/featureCountsTables/EcC3000_5_min_featureCounts_TPM_normalized.tsv",
       sep = "\t", quote = F, row.names = F)
