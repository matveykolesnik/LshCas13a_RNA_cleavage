library(Gviz)
library(rtracklayer)
options(ucscChromosomeNames = F)

setwd("~/data/LshCas13a_RNA_cleavage/LshCas13a_M13_infection")

M13_SeqID <- "M13+Sp1rev_C-PFS"
AnnotationFile <- "Annotations/M13_annotation.gff"
Annotation.GR <- readGFFAsGRanges(AnnotationFile)

### load coverage table
library(tidyverse)

ReadsCoverageTable <- "Results/Tables/MergedCoverageCounts/Reads_coverage_merged.tsv.gz"
ReadsCoverageTable.df <- read.delim(ReadsCoverageTable, sep = "\t", stringsAsFactors = F) %>% 
  filter(SeqID == M13_SeqID)

ReadsCoverageTable_forward.m <- as.matrix(ReadsCoverageTable.df[ReadsCoverageTable.df$Strand == "+",][4:ncol(ReadsCoverageTable.df)])
ReadsCoverageTable_reverse.m <- as.matrix(ReadsCoverageTable.df[ReadsCoverageTable.df$Strand == "-",][4:ncol(ReadsCoverageTable.df)])

reads_per_lib <- "Results/Tables/reads_per_lib.tsv"
reads_per_lib.df <- read.delim(reads_per_lib, sep="\t", header = F, col.names = c("libname", "size")) %>% 
  mutate(libname = gsub("_sorted.bam", "", libname),
         size = size/10**6) %>% 
  arrange(factor(libname, levels = colnames(ReadsCoverageTable.df[4:ncol(ReadsCoverageTable.df)])))

#normalization
ReadsCoverageTable_forward.m.norm <- sweep(ReadsCoverageTable_forward.m, 2, reads_per_lib.df$size, "/")
ReadsCoverageTable_forward.m.norm[ReadsCoverageTable_forward.m.norm == 0] <- 0.1
ReadsCoverageTable_forward.m.norm <- log10(ReadsCoverageTable_forward.m.norm)

ReadsCoverageTable_reverse.m.norm <- sweep(ReadsCoverageTable_reverse.m, 2, reads_per_lib.df$size, "/")
ReadsCoverageTable_reverse.m.norm[ReadsCoverageTable_reverse.m.norm == 0] <- 0.1
ReadsCoverageTable_reverse.m.norm <- log10(ReadsCoverageTable_reverse.m.norm)

max_value <- max(c(ReadsCoverageTable_forward.m.norm, ReadsCoverageTable_reverse.m.norm))
min_value <- min(c(ReadsCoverageTable_forward.m.norm, ReadsCoverageTable_reverse.m.norm))

#oh blyat
T_samples <- c("M13_T1_S7", "M13_T2_S8", "M13_T3_S9")
NT_samples <- c("M13_NT1_S10", "M13_NT2_S11", "M13_NT3_S12")

M13_annotation_track <- AnnotationTrack(start = Annotation.GR@ranges@start, 
                                        width = Annotation.GR@ranges@width,
                                        chromosome = M13_SeqID,
                                        strand = as.character(Annotation.GR@strand),
                                        id = Annotation.GR$ID,
                                        stacking="dense",
                                        fill = "grey",
                                        col = "black",
                                        lex = 1,
                                        lty = 1,
                                        shape = "fixedArrow",
                                        arrowHeadWidth = 10, 
                                        lwd = 1,
                                        fontcolor.feature="black",
                                        rotation.item = 90,
                                        cex = 0.5,
                                        labelPos = "below")

axisTrack <- GenomeAxisTrack(cex = 1, labelPos="below", size = 1)

ReadsCoverageTable_forward_T.dTrack <-  DataTrack(start = seq(1, nrow(ReadsCoverageTable_forward.m.norm), 1),
                                                  width = 1,
                                                  chromosome = M13_SeqID,
                                                  genome = "M13",
                                                  name = "Targeting samples,\nlg(coverage)",
                                                  data = t(ReadsCoverageTable_forward.m.norm[,T_samples]),
                                                  col = "red",
                                                  ylim = c(min_value, max_value))
ReadsCoverageTable_forward_NT.dTrack <- DataTrack(start = seq(1, nrow(ReadsCoverageTable_forward.m.norm), 1),
                                                  width = 1,
                                                  chromosome = M13_SeqID,
                                                  genome = "M13",
                                                  name = "Nontargeting samples,\nlg(coverage)",
                                                  data = t(ReadsCoverageTable_forward.m.norm[,NT_samples]),
                                                  col = "grey",
                                                  ylim = c(min_value, max_value))

# ReadsCoverageTable_reverse_T.dTrack <- DataTrack(start = seq(1, nrow(ReadsCoverageTable_reverse.m.norm), 1),
#                                                  width = 1,
#                                                  chromosome = M13_SeqID,
#                                                  genome = "M13",
#                                                  name = "Targeting samples,\ncoverage",
#                                                  data = t(ReadsCoverageTable_reverse.m.norm[,T_samples]),
#                                                  col = "red",
#                                                  ylim = c(0, max_value))
# ReadsCoverageTable_reverse_NT.dTrack <- DataTrack(start = seq(1, nrow(ReadsCoverageTable_reverse.m.norm), 1),
#                                                   width = 1,
#                                                   chromosome = M13_SeqID,
#                                                   genome = "M13",
#                                                   name = "Nontargeting samples,\ncoverage",
#                                                   data = t(ReadsCoverageTable_reverse.m.norm[,NT_samples]),
#                                                   col = "red",
#                                                   ylim = c(0, max_value))

svg(filename = "Results/Pictures/M13_reads_coverage_lg_transformed.svg", width = 15, height = 7)
plotTracks(c(ReadsCoverageTable_forward_T.dTrack, ReadsCoverageTable_forward_NT.dTrack, 
             M13_annotation_track, 
#             ReadsCoverageTable_reverse_T.dTrack, ReadsCoverageTable_reverse_NT.dTrack,
             axisTrack),
           featureAnnotation = "id", type = c("a", "confint"))
dev.off()

