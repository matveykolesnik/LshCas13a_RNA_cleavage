library(dplyr)
library(tidyverse)
library(edgeR)
library(ggfittext)
library(reshape2)
library(Biostrings)
library(data.table)

RefSeqs <- readDNAStringSet("Reference_sequences/Nontargeting_samples_refseqs.fasta", format = "fasta")
chromosome = RefSeqs$NC_000913.3
pC008 = RefSeqs$pC008

targeting_samples <- c("IJ01_S1", "IJ02_S2", "IJ03_S3")
nontargeting_samples <- c("IJ04_S4", "IJ05_S5", "IJ06_S6")

N5E_tables <- read.delim("Results/Tables/Merged_ends_counts/N5E_TAP_untreated_T_vs_NT.tsv.gz", sep = "\t")

N5E_tables$T_aveLogCPM <- aveLogCPM(N5E_tables[targeting_samples], prior.count = 1)
N5E_tables$NT_aveLogCPM <- aveLogCPM(N5E_tables[nontargeting_samples], prior.count = 1)
N5E_tables <- N5E_tables %>% 
  mutate(., aveLogFC = T_aveLogCPM - NT_aveLogCPM)

N5E_tables[4:9] <- cpm(N5E_tables[4:9])

fwrite(x = N5E_tables, file = "Results/Tables/Ends_counts_normalized/N5E_CPM.tsv", sep = "\t", quote = F, row.names = F)
# write.table(N5E_tables, "Results/Tables/Ends_counts_normalized/N5E_CPM.tsv", sep = "\t", quote = F, row.names = F)

#rpmH TCS
rpmH_TCS <- 3884341
rpmH_start <- 3884336

rpmH_TCS_surroundings <- N5E_tables[N5E_tables$SeqID == "NC_000913.3" & N5E_tables$Strand == "+" & N5E_tables$Pos >= rpmH_TCS-5 & N5E_tables$Pos <= rpmH_TCS+15, ] %>% 
  rowwise() %>% 
  mutate(., Letter = as.character(chromosome[Pos])) %>% 
  mutate(., Pos = Pos-rpmH_start+1) %>% 
  mutate(., Letter = sub("T", "U", Letter))

write.table(rpmH_TCS_surroundings, "Results/Tables/Picture_source_tables/Ecd10LVM_rpmH_TCS_raw.tsv", sep = "\t", quote = F, row.names = F)

rpmH_TCS_surroundings.melted <- reshape2::melt(rpmH_TCS_surroundings, id.vars=c("SeqID", "Pos", "Strand", "Letter"), variable.name = "Exp_group", value.name = "N5E_CPM") %>% 
  mutate(., Exp_group = ifelse(Exp_group %in% targeting_samples, "Targeting", "Nontargeting"))

rpmH_TCS_surroundings.melted.agg <- aggregate(N5E_CPM~SeqID*Pos*Strand*Letter*Exp_group, data=rpmH_TCS_surroundings.melted, FUN="mean")
rpmH_TCS_surroundings.melted.agg$SD <- aggregate(N5E_CPM~SeqID*Pos*Strand*Letter*Exp_group, data=rpmH_TCS_surroundings.melted, FUN="sd")$N5E_CPM

write.table(rpmH_TCS_surroundings.melted.agg, "Results/Tables/Picture_source_tables/Ecd10LVM_rpmH_TCS.tsv", sep = "\t", quote = F, row.names = F)

#plot coverage in surroundings of bla TCS
bla_TCS <- 1576
bla_start = 1565

bla_TCS_surroundings <- N5E_tables[N5E_tables$SeqID == "pC008" & N5E_tables$Strand == "+" & N5E_tables$Pos >= bla_TCS-11 & N5E_tables$Pos <= bla_TCS+9, ] %>% 
  rowwise() %>% 
  mutate(., Letter = as.character(pC008[Pos])) %>% 
  mutate(., Pos = Pos-bla_start+1) %>% 
  mutate(., Letter = sub("T", "U", Letter))

write.table(bla_TCS_surroundings, "Results/Tables/Picture_source_tables/Ecd10LVM_bla_TCS_raw.tsv", sep = "\t", quote = F, row.names = F)

bla_TSS <- 1530

bla_TCS_surroundings_extended <- N5E_tables[N5E_tables$SeqID == "pC008" & N5E_tables$Strand == "+" & N5E_tables$Pos >= bla_TSS & N5E_tables$Pos <= bla_TCS+9, ] %>% 
  rowwise() %>% 
  mutate(., Letter = as.character(pC008[Pos])) %>% 
  mutate(., Pos = Pos-bla_start+1) %>% 
  mutate(., Letter = sub("T", "U", Letter))

write.table(bla_TCS_surroundings_extended, "Results/Tables/Picture_source_tables/Ecd10LVM_bla_TCS_surroundings_raw.tsv", sep = "\t", quote = F, row.names = F)

####
bla_TCS_surroundings.melted <- reshape2::melt(bla_TCS_surroundings, id.vars=c("SeqID", "Pos", "Strand", "Letter"), variable.name = "Exp_group", value.name = "N5E_CPM") %>% 
  mutate(., Exp_group = ifelse(Exp_group %in% targeting_samples, "Targeting", "Nontargeting"))

bla_TCS_surroundings.melted.agg <- aggregate(N5E_CPM~SeqID*Pos*Strand*Letter*Exp_group, data=bla_TCS_surroundings.melted, FUN="mean")
bla_TCS_surroundings.melted.agg$SD <- aggregate(N5E_CPM~SeqID*Pos*Strand*Letter*Exp_group, data=bla_TCS_surroundings.melted, FUN="sd")$N5E_CPM

write.table(bla_TCS_surroundings.melted.agg, "Results/Tables/Picture_source_tables/Ecd10LVM_bla_TCS.tsv", sep = "\t", quote = F, row.names = F)

# bla_TCS_plot <- ggplot(bla_TCS_surroundings.melted.agg, aes(x=Pos, y=N5E_CPM)) +
#   geom_bar(stat = "identity", fill="grey") +
#   geom_errorbar(aes(ymin=N5E_CPM-SD, ymax=N5E_CPM+SD), position="dodge", width=0.5) +
#   geom_point(data = bla_TCS_surroundings.melted, aes(Pos, N5E_CPM)) +
#   xlab("Nucleotide position of bla gene") +
#   ylab("Mean value of CPM-normalized N5E") +
#   ggtitle("N5E coverage across the beginning of bla gene") +
#   facet_wrap(~Exp_group, nrow = 2) +
#   theme_classic() +
#   theme(axis.title = element_text(size = 14),
#         axis.text = element_text(size = 14),
#         title = element_text(size = 14),
#         strip.text = element_text(size = 14))
# 
# ggsave("Results/Pictures/WT_bla_TCS.png", plot = bla_TCS_plot, device = "png", scale = 1.1)
# ggsave("Results/Pictures/WT_bla_TCS.svg", plot = bla_TCS_plot, device = "svg", scale = 1.1)

#plot N5E coverage across lysT tRNA gene
lysT_start = 780554
lysT_end = 780629
lysT_AC_start = 780587-lysT_start

lysT_surroundings <- N5E_tables[N5E_tables$SeqID == "NC_000913.3" & N5E_tables$Strand == "+" & N5E_tables$Pos >= lysT_start & N5E_tables$Pos <= lysT_end, ] %>% 
  rowwise() %>% 
  mutate(., Letter = as.character(chromosome[Pos])) %>% 
  mutate(., Pos = Pos-lysT_start+1) %>% 
  mutate(., Letter = sub("T", "U", Letter))

write.table(lysT_surroundings, "Results/Tables/Picture_source_tables/Ecd10LVM_lysT_TCS_raw.tsv", sep = "\t", quote = F, row.names = F)

lysT_surroundings.melted <- reshape2::melt(lysT_surroundings, id.vars=c("SeqID", "Pos", "Strand", "Letter"), variable.name = "Exp_group", value.name = "N5E_CPM") %>% 
  mutate(., Exp_group = ifelse(Exp_group %in% targeting_samples, "Targeting", "Nontargeting"))

lysT_surroundings.melted.agg <- aggregate(N5E_CPM~SeqID*Pos*Strand*Letter*Exp_group, data=lysT_surroundings.melted, FUN="mean")
lysT_surroundings.melted.agg$SD <- aggregate(N5E_CPM~SeqID*Pos*Strand*Letter*Exp_group, data=lysT_surroundings.melted, FUN="sd")$N5E_CPM

lysT_surroundings.melted_var.agg <- aggregate(N5E_CPM~SeqID*Pos*Strand*Letter*Exp_group, data=lysT_surroundings.melted, FUN="mean")
lysT_surroundings.melted_var.agg$VAR <- aggregate(N5E_CPM~SeqID*Pos*Strand*Letter*Exp_group, data=lysT_surroundings.melted, FUN="var")$N5E_CPM

write.table(lysT_surroundings.melted.agg, 'Results/Tables/Picture_source_tables/Ecd10LVM_lysT_TCS.tsv', sep = "\t", quote = F, row.names = F)
#write.table(lysT_surroundings.melted_var.agg, "../Drafts/Source_tables/WT_lysT_TCS_mean_var.tsv", sep = "\t", quote = F, row.names = F)

# lysT_TCS_plot <- ggplot(lysT_surroundings.melted.agg, aes(x=Pos, y=N5E_CPM)) +
#   geom_bar(stat = "identity", fill="grey") +
#   geom_errorbar(aes(ymin=N5E_CPM-SD, ymax=N5E_CPM+SD), position="dodge", width=0.5) +
#   geom_point(data = lysT_surroundings.melted, aes(Pos, N5E_CPM)) +
#   xlab("Nucleotide position of lysT gene") +
#   ylab("Mean value of CPM-normalized N5E") +
#   ggtitle("N5E coverage across lysT tRNA gene") +
#   facet_wrap(~Exp_group, nrow = 2) +
#   theme_classic() +
#   geom_rect(aes(xmin=lysT_AC_start-0.5, ymin=-40, xmax=lysT_AC_start+2.5, ymax=-10), fill="lightblue", colour="black") +
#   geom_fit_text(aes(xmin=lysT_AC_start-0.5, ymin=-40, xmax=lysT_AC_start+2.5, ymax=-10, label="AC")) +
#   theme(axis.title = element_text(size = 14),
#         axis.text = element_text(size = 14),
#         title = element_text(size = 14),
#         strip.text = element_text(size = 14))
# 
# ggsave("Results/Pictures/WT_lysT_TCS.png", plot = lysT_TCS_plot, device = "png", scale = 1.5)
# ggsave("Results/Pictures/WT_lysT_TCS.svg", plot = lysT_TCS_plot, device = "svg", scale = 1.5)

#plot rrsH TCS
rrsH_TCS = 225267
rrsH_start = 223771
rrsH_end = 225312

anti_SD_start = 225305-rrsH_start
anti_SD_end = 225310-rrsH_start

rrsH_TCS_surroundings <- N5E_tables[N5E_tables$SeqID == "NC_000913.3" & N5E_tables$Strand == "+" & N5E_tables$Pos >= rrsH_TCS-10 & N5E_tables$Pos <= rrsH_end, ] %>% 
  rowwise() %>% 
  mutate(., Letter = as.character(chromosome[Pos])) %>% 
  mutate(., Pos = Pos-rrsH_start) %>% 
  mutate(., Letter = sub("T", "U", Letter))

write.table(rrsH_TCS_surroundings, "Results/Tables/Picture_source_tables/Ecd10LVM_rrsH_TCS_raw.tsv", sep = "\t", quote = F, row.names = F)

rrsH_TCS_surroundings.melted <- reshape2::melt(rrsH_TCS_surroundings, id.vars=c("SeqID", "Pos", "Strand", "Letter"), variable.name = "Exp_group", value.name = "N5E_CPM") %>% 
  mutate(., Exp_group = ifelse(Exp_group %in% targeting_samples, "Targeting", "Nontargeting"))

rrsH_TCS_surroundings.melted.agg <- aggregate(N5E_CPM~SeqID*Pos*Strand*Letter*Exp_group, data=rrsH_TCS_surroundings.melted, FUN="mean")
rrsH_TCS_surroundings.melted.agg$SD <- aggregate(N5E_CPM~SeqID*Pos*Strand*Letter*Exp_group, data=rrsH_TCS_surroundings.melted, FUN="sd")$N5E_CPM

write.table(rrsH_TCS_surroundings.melted.agg, 'Results/Tables/Picture_source_tables/Ecd10LVM_rrsH_TCS.tsv', sep = "\t", quote = F, row.names = F)

# rrsH_TCS_plot <- ggplot(rrsH_TCS_surroundings.melted.agg, aes(x=Pos, y=N5E_CPM)) +
#   geom_bar(stat = "identity", fill="grey") +
#   geom_errorbar(aes(ymin=N5E_CPM-SD, ymax=N5E_CPM+SD), position="dodge", width=0.5) +
#   geom_point(data = rrsH_TCS_surroundings.melted, aes(Pos, N5E_CPM)) +
#   xlab("Nucleotide position of rrsH gene") +
#   ylab("Mean value of CPM-normalized N5E") +
#   ggtitle("N5E coverage across the end of rrsH gene") +
#   facet_wrap(~Exp_group, nrow = 2) +
#   theme_classic() +
#   geom_rect(aes(xmin=anti_SD_start-0.5, ymin=-2, xmax=anti_SD_end+0.5, ymax=-0.5), fill="lightblue", colour="black") +
#   geom_fit_text(aes(xmin=anti_SD_start-0.5, ymin=-2, xmax=anti_SD_end+0.5, ymax=-0.5, label="anti-SD")) +
#   theme(axis.title = element_text(size = 14),
#         axis.text = element_text(size = 14),
#         title = element_text(size = 14),
#         strip.text = element_text(size = 14))
# 
# ggsave("Results/Pictures/WT_rrsH_TCS.png", plot = rrsH_TCS_plot, device = "png")
# ggsave("Results/Pictures/WT_rrsH_TCS.svg", plot = rrsH_TCS_plot, device = "svg")

#save(bla_TCS_plot, lysT_TCS_plot, rrsH_TCS_plot, file = "Results/Pictures/WT_TCS_examples.R")

