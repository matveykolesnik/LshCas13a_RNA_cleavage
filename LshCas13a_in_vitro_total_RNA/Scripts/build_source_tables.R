library(dplyr)
library(tidyr)
library(Biostrings)

setwd("~/data/LshCas13a_RNA_cleavage/LshCas13a_in_vitro_total_RNA/")

#build source tables for tRNAs
RefSeqs <- readDNAStringSet("Reference_sequences/NC_000913.3.fasta", format = "fasta")
chromosome = RefSeqs$NC_000913.3

tRNA_annotation <- read.delim("Annotations/NC_000913.3_Aragorn_tRNAs_renamed.tsv", 
                              sep = "\t", col.names = c("tRNA_type", "tRNA_start", "tRNA_end", "tRNA_strand", "AC_pos", "AC"))

tRNA_tables_dir <- "Results/Tables/Source_tables/tRNA_source_tables"
dir.create(tRNA_tables_dir, recursive = T)

build_tRNA_table <- function(N5E_table, tRNA_type, tRNA_start, tRNA_end, tRNA_strand, AC_pos, output_dir) {
  print(tRNA_type)
  tRNA_table <- N5E_table[N5E_table$SeqID == "NC_000913.3" & N5E_table$Strand == tRNA_strand & N5E_table$Pos >= tRNA_start & N5E_table$Pos <= tRNA_end, ]
  
  if (tRNA_strand == "+") {
    tRNA_table <- tRNA_table %>% 
      rowwise() %>% 
      mutate(., Letter = as.character(chromosome[Pos])) %>% 
      mutate(., Letter = sub("T", "U", Letter)) %>% 
      mutate(., Pos = Pos-tRNA_start+1) %>% 
      arrange(., Pos)
    
  } else {
    tRNA_table <- tRNA_table %>% 
      rowwise() %>% 
      mutate(., Letter = as.character(complement(chromosome)[Pos])) %>% 
      mutate(., Letter = sub("T", "U", Letter)) %>% 
      mutate(., Pos = tRNA_end-Pos+1) %>% 
      arrange(., Pos)
  }
  AC_coords <- seq(AC_pos, AC_pos+2)
  
  output_file <- file.path(output_dir, paste0(tRNA_type, "_", paste(c(tRNA_start, tRNA_end, AC_coords), sep = "_", collapse = "_"), ".tsv"))
  write.table(tRNA_table, file = output_file, sep = "\t", quote = F, row.names = F)
}

N5E_counts_table <- read.delim("Results/Tables/Merged_ends_counts/N5E_T_vs_NT_CPM_normalized.tsv.gz",
                               sep = "\t", 
                               header = T, 
                               stringsAsFactors = F)
apply(tRNA_annotation, 1, function(x) build_tRNA_table(N5E_table = N5E_counts_table,
                                                       tRNA_type = x["tRNA_type"], 
                                                       tRNA_start = as.integer(x["tRNA_start"]), 
                                                       tRNA_end = as.integer(x["tRNA_end"]), 
                                                       tRNA_strand = x["tRNA_strand"], 
                                                       AC_pos = as.integer(x["AC_pos"]),
                                                       output_dir = tRNA_tables_dir))

#create source table for rrsH transcript