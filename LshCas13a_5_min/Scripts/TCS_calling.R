library(dplyr)
library(data.table)
library(tidyr)
library(edgeR)

setwd("/home/ledum/data/LshCas13a_RNA_cleavage/LshCas13a_5_min")
dir.create("Results/Tables/TCS_detection_tables")

N5E_counts_file <- "Results/Tables/Merged_ends_counts/N5E_T_vs_NT.tsv.gz"
design_file <- "design.tsv"

#load table
N5E_counts_table <- read.delim(N5E_counts_file, sep = "\t", stringsAsFactors = F, header = T)
#rename rows; each rowname looks like "SeqID:Pos:Strand"
rownames(N5E_counts_table) <- apply(N5E_counts_table, 1, function(x) paste0(x["SeqID"], ":", x["Pos"], ":", x["Strand"]))

design_table <- read.delim(design_file, sep = "\t", stringsAsFactors = F, header = T)
samples <- design_table[match(colnames(N5E_counts_table)[4:ncol(N5E_counts_table)], design_table$Sample),]$Exp_group
design <- model.matrix(~samples)

#convert data frame to DGE list
dge <- DGEList(counts = N5E_counts_table[4:9], group = samples)
#dge <- dge[rowSums(cpm(dge)>=1) >= 3, , keep.lib.sizes=FALSE]
dge <- dge[rowSums(cpm(dge)[,design_table[design_table$Exp_group == "targeting", "Sample"]]>=10) >= 3, , keep.lib.sizes=FALSE]
#dge <- dge[rowSums(cpm(dge)>=1) == 6, , keep.lib.sizes=FALSE]
#calculate factors of normalization
dge <- calcNormFactors(dge, method = "TMM")
#estimate dispersion
dge <- estimateDisp(dge, design)

plotMDS(dge, main="edgeR MDS Plot", col=c(rep("blue", 3), rep("red", 3)), pch = 19,
        xlab = "Dimension 1",
        ylab = "Dimension 2")
legend("topleft", inset=.02, c("Targeting samples","Non-targeting samples"), fill=c("red", "blue"), horiz=FALSE, cex=0.8)

fit <- glmFit(dge, design = design)
lrt <- glmLRT(fit, coef = 2)

result_table <- as.data.table(lrt$table, keep.rownames = "feature") %>% 
  separate(., col = "feature", into = c("SeqID", "Pos", "Strand"), sep = ":") %>% 
  mutate(., Pos = as.integer(Pos),
         PValue.adj = p.adjust(PValue, method = "BH"))

write.table(result_table, file = "Results/Tables/TCS_detection_tables/LRTest_table.tsv", sep = "\t", row.names = F, quote = F)

library(rtracklayer)
result_table.GR <- with(result_table[result_table$PValue.adj <= 10**-6 & result_table$logFC>0, ], GRanges(seqnames = SeqID, strand = Strand, ranges = IRanges(Pos, Pos), ID = paste(SeqID, Pos, logFC, PValue.adj, sep = ":")))

export(result_table.GR, "Annotations/predicted_TCS.gff3")
