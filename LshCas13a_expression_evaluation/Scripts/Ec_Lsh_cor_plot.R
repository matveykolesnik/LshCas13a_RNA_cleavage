library(ggplot2)
library(ggrepel)
library(data.table)
library(dplyr)

Ec_featureCounts.TPM <- "Results/Tables/featureCountsTables/Ec_Shmakov_data_TPM_normalized.tsv"
Lsh_featureCounts.TPM <- "Results/Tables/featureCountsTables/Lsh_Shmakov_data_TPM_normalized.tsv"

Ec_featureCounts.TPM.dt <- fread(Ec_featureCounts.TPM) %>% 
  select(Name, EcTPM = E_coli_data.bam)

Lsh_featureCounts.TPM.dt <- fread(Lsh_featureCounts.TPM) %>% 
  select(Name, LshTPM = L_shahii_data.bam) %>% 
  mutate(Name = ifelse(Name == "cas13a", "LshCas13a", Name))

Ec_Lsh.TPM <- inner_join(Ec_featureCounts.TPM.dt, Lsh_featureCounts.TPM.dt, by="Name") %>% 
  filter_at(c("EcTPM", "LshTPM"), all_vars(. > 0))
print(nrow(Ec_Lsh.TPM))

labeled_genes <- c("LshCas13a", "ssrA", "ffs", "rnpB", "ssrS")

Ec_Lsh.TPM.labeled <- Ec_Lsh.TPM %>% 
  mutate(ref = ifelse(Name %in% labeled_genes, Name, "")) %>% 
  mutate(ref = ifelse(Name == "LshSpacer_1", "crRNA", ref))

cor.test(x = log10(Ec_Lsh.TPM.labeled$EcTPM), y = log10(Ec_Lsh.TPM.labeled$LshTPM))

correlation_plot <- ggplot(Ec_Lsh.TPM.labeled, aes(x = log10(EcTPM), y = log10(LshTPM))) +
  geom_point()+
  geom_smooth(method=lm) +
  #scale_color_manual(values = c("black", "red")) +
  geom_text_repel(aes(label = ref)) +
  xlab("logTPM of E. coli gene transcripts") +
  ylab("logTPM of L. shahii gene transcripts") +
  theme_bw()

ggsave(filename = "Results/Pictures/Eco_Lsh_cor_plot_cas13a_added.png", plot = correlation_plot, dpi = "retina", height = 10, width = 10)

