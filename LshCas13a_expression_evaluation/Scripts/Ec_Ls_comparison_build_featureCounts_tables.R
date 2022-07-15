library(Rsubread)

EcBAMFile <- "Alignments/Shmakov_et_al/E_coli_RNASeq/E_coli_data.bam"
LshBAMfile <- "Alignments/Shmakov_et_al/L_shahii_RNASeq/L_shahii_data.bam"

EcAnnotation <- "Annotations/E_coli_merged_annotation_spacer.gff3"
LshAnnotation <- "Annotations/L_shahii_DSM_19757_spacer.gff3"

Ec_featureCounts <- featureCounts(files = EcBAMFile,
                                  annot.ext = EcAnnotation,
                                  isGTFAnnotationFile = T,
                                  GTF.featureType = "gene",
                                  GTF.attrType = "ID",
                                  GTF.attrType.extra = "Name",
                                  allowMultiOverlap = T,
                                  strandSpecific = 1,
                                  isPairedEnd = T,
                                  countReadPairs = T,
                                  requireBothEndsMapped = F,
                                  nthreads = 16)

Lsh_featureCounts <- featureCounts(files = LshBAMfile,
                                   annot.ext = LshAnnotation,
                                   isGTFAnnotationFile = T,
                                   GTF.featureType = "gene",
                                   GTF.attrType = "ID",
                                   GTF.attrType.extra = "Name",
                                   allowMultiOverlap = T,
                                   strandSpecific = 1,
                                   isPairedEnd = T,
                                   countReadPairs = T,
                                   requireBothEndsMapped = F,
                                   nthreads = 16)

library(data.table)
library(dplyr)

Ec_featureCounts_raw <-  left_join(as.data.table(Ec_featureCounts$counts, keep.rownames = "GeneID"),
                                   Ec_featureCounts$annotation,
                                   by="GeneID") %>%
  select(-c("Start", "End", "Strand", "Length", "Chr")) %>% 
  select(GeneID, Name, everything())
fwrite(x = Ec_featureCounts_raw, file = "Results/Tables/featureCountsTables/Ec_Shmakov_data_raw_counts.tsv",
       sep = "\t", quote = F, row.names = F)

Lsh_featureCounts_raw <-  left_join(as.data.table(Lsh_featureCounts$counts, keep.rownames = "GeneID"),
                                    Lsh_featureCounts$annotation,
                                    by="GeneID") %>%
  select(-c("Start", "End", "Strand", "Length", "Chr")) %>% 
  select(GeneID, Name, everything())
fwrite(x = Lsh_featureCounts_raw, file = "Results/Tables/featureCountsTables/Lsh_Shmakov_data_raw_counts.tsv",
       sep = "\t", quote = F, row.names = F)

Ec_featureCounts_matrix <- as.matrix(Ec_featureCounts$counts)
Ec_featureCounts_matrix.norm_by_length <- Ec_featureCounts_matrix/Ec_featureCounts$annotation$Length
Ec_featureCounts.TPM <- t(t(Ec_featureCounts_matrix.norm_by_length)*10**6/colSums(Ec_featureCounts_matrix.norm_by_length)) %>% 
  as.data.table(., keep.rownames = "GeneID") %>% 
  left_join(., Ec_featureCounts$annotation, by="GeneID") %>% 
  select(-c("Start", "End", "Strand", "Length", "Chr")) %>% 
  select(GeneID, Name, everything())

fwrite(x = Ec_featureCounts.TPM, file = "Results/Tables/featureCountsTables/Ec_Shmakov_data_TPM_normalized.tsv",
       sep = "\t", quote = F, row.names = F)

Lsh_featureCounts_matrix <- as.matrix(Lsh_featureCounts$counts)
Lsh_featureCounts_matrix.norm_by_length <- Lsh_featureCounts_matrix/Lsh_featureCounts$annotation$Length
Lsh_featureCounts.TPM <- t(t(Lsh_featureCounts_matrix.norm_by_length)*10**6/colSums(Lsh_featureCounts_matrix.norm_by_length)) %>% 
  as.data.table(., keep.rownames = "GeneID") %>% 
  left_join(., Lsh_featureCounts$annotation, by="GeneID") %>% 
  select(-c("Start", "End", "Strand", "Length", "Chr")) %>% 
  select(GeneID, Name, everything())

fwrite(x = Lsh_featureCounts.TPM, file = "Results/Tables/featureCountsTables/Lsh_Shmakov_data_TPM_normalized.tsv",
       sep = "\t", quote = F, row.names = F)
