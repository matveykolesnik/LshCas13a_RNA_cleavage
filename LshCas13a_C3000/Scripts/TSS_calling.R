library(dplyr)
library(data.table)
library(tidyr)
library(edgeR)

setwd("~/data/LshCas13a_RNA_cleavage/LshCas13a_C3000/")
dir.create("Results/Tables/TSS_detection_tables")

N5E_counts_file <- "Results/Tables/Merged_ends_counts/N5E_T_TAP_treated_vs_TAP_untreated.tsv.gz"
design_file <- "design.tsv"

N5E_counts_table <- read.delim(N5E_counts_file, sep = "\t", header = T, stringsAsFactors = F)
rownames(N5E_counts_table) <- apply(N5E_counts_table, 1, function(x) paste0(x["SeqID"], ":", x["Pos"], ":", x["Strand"]))

design_table <- read.delim(design_file, sep = "\t", stringsAsFactors = F, header = T)
samples <- factor(design_table[match(colnames(N5E_counts_table)[4:ncol(N5E_counts_table)], design_table$Sample),]$TAP_treatment)
design <- model.matrix(~0+samples)
colnames(design) <- levels(samples)

TAPvsNoTAP <- makeContrasts(treated-untreated, levels = design)

dge <- DGEList(counts = N5E_counts_table[4:9], group = samples)
dge <- dge[rowSums(cpm(dge)>=10) >= 3, , keep.lib.sizes=FALSE]
dge <- calcNormFactors(dge, method = "TMM")
dge <- estimateDisp(dge, design)

plotMDS(dge, main="edgeR MDS Plot", col=c(rep("blue", 3), rep("red", 3)), pch = 19,
        xlab = "Dimension 1",
        ylab = "Dimension 2")
legend("topleft", inset=.02, c("Targeting samples","Non-targeting samples"), fill=c("red", "blue"), horiz=FALSE, cex=0.8)

plotBCV(dge)

fit <- glmFit(dge, design = design)
lrt <- glmLRT(fit, contrast = TAPvsNoTAP)

result_table <- as.data.table(lrt$table, keep.rownames = "feature") %>% 
  separate(., col = "feature", into = c("SeqID", "Pos", "Strand"), sep = ":") %>% 
  mutate(., Pos = as.integer(Pos),
         PValue.adj = p.adjust(PValue, method = "BH")) %>% 
  arrange(., PValue.adj)

write.table(result_table, file = "Results/Tables/TSS_detection_tables/T_cells_TSS_predictions.tsv", sep = "\t", row.names = F, quote = F)
