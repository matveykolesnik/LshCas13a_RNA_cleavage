library(dplyr)
library(data.table)
library(rtracklayer)

setwd("/home/ledum/data/LshCas13a_RNA_cleavage/LshCas13a_C3000/")

annotation_file <- "Annotations/Summarized_annotation.gff3"

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

overlaps_table.merged <- inner_join(x = TCS_df, y = overlaps_table, by = c("SeqID", "Pos", "Strand")) %>% 
  arrange(., PValue.adj)

write.table(overlaps_table.merged, "Results/Tables/TCS_detection_tables/LRTest_table_and_genome_features.tsv", sep = "\t", row.names = F, quote = F)

library(ggplot2)
library(ggrepel)

overlaps_table.merged.viz <- overlaps_table.merged
overlaps_table.merged.viz$feature_label <- as.character("")
t = 10
overlaps_table.merged.viz[1:t,]$feature_label <- overlaps_table.merged.viz[1:t,]$MatchedFeatureGene

volcano_plot <- ggplot(overlaps_table.merged.viz, aes(x=logFC, y=-log10(PValue.adj), colour=MatchedFeatureType)) +
  geom_point() +
  geom_text_repel(aes(label=feature_label)) +
  ggtitle(label = "EcC3000 + M13 TCS volcano plot") +
  xlab(label = "log2 fold change between targeting and nontargeting samples") +
  ylab(label = "-log10 adjusted p-value") +
  guides(color=guide_legend("Type of genomic feature")) +
  theme_bw()

ggsave("Results/Pictures/EcC3000_M13_volcano_draft.png", plot = volcano_plot, dpi = 300)
