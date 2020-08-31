library(dplyr)
library(data.table)
library(rtracklayer)

setwd("/home/matvey/data/LshCas13a_RNA_cleavage/LshCas13a_in_vitro_tRNAs/")

annotation_file <- "Annotations/NC_000913.3.gff3"

#loading annotation
#annotation_GR <- readGFFAsGRanges(filepath = annotation_file, filter = list("type"=c("CDS", "rRNA", "tRNA")))
annotation_GR <- readGFFAsGRanges(filepath = annotation_file)

#loading TCS table
TCS_df <- read.delim(file = "Results/Tables/TCS_detection_tables/LRTest_table.tsv", sep = "\t", stringsAsFactors = F)
#generating GRanges object from TCS table
TCS_df.GR <- with(TCS_df, GRanges(seqnames = SeqID, strand = Strand, ranges = IRanges(Pos, Pos)))
#finding overlaps between genomic features and TCS
overlaps <- findOverlaps(query = TCS_df.GR, subject = annotation_GR, type = "within")
#creating table of overlaps

overlaps_table <- data.table(SeqID = as.character(TCS_df.GR[overlaps@from]@seqnames),
                             Pos = as.integer(TCS_df.GR[overlaps@from]@ranges@start),
                             Strand = as.character(TCS_df.GR[overlaps@from]@strand),
                             MatchedFeatureType = annotation_GR[overlaps@to]$type,
                             MatchedFeatureID = annotation_GR[overlaps@to]$ID,
                             MatchedFeatureGene = annotation_GR[overlaps@to]$gene,
                             MatchedFeatureDescription = annotation_GR[overlaps@to]$product,
                             MatchedFeatureStart = annotation_GR[overlaps@to]@ranges@start,
                             MatchedFeatureEnd = annotation_GR[overlaps@to]@ranges@start + annotation_GR[overlaps@to]@ranges@width-1,
                             MatchedFeatureStrand = as.character(annotation_GR[overlaps@to]@strand)) %>% 
  filter(., !(MatchedFeatureType %in% c("gene", "exon")))

overlaps_table.merged <- inner_join(x = TCS_df, y = overlaps_table, by = c("SeqID", "Pos", "Strand"))

write.table(overlaps_table.merged, "Results/Tables/TCS_detection_tables/LRTest_table_and_genome_features.tsv", sep = "\t", row.names = F, quote = F)