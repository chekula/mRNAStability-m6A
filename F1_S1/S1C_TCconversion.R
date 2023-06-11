library (tidyverse)
library (dplyr)
library (ggplot2)
library(gridExtra)
library(ggrepel)
library(rstatix)


## load conversion rates for the whole cell

ftc <- read.csv('../data/ftc_raw_whole_cell.csv') %>% 
  select(gene_id, condition, replicate, timepoint, ftc)

pseudoR2 <- read.csv('../data/TtoC_full_results.csv') %>% 
  select(gene_id, condition, replicate, pseudoR2)

conversions <- ftc %>% left_join(pseudoR2, by=c('gene_id', 'condition', 'replicate')) %>%
  filter(pseudoR2>0.95) %>% 
  select(-pseudoR2) %>% 
group_by(gene_id, condition, timepoint) %>%
  mutate(ftc_mean = mean(ftc, na.rm = TRUE)) %>%
  select(-c(replicate, ftc)) %>%
  filter(ftc_mean != 0) %>% 
  distinct() %>%
  ungroup()

gfp_caf <- conversions %>%
  #filter(condition == "GFP") %>% 
  mutate(timepoint = as.factor(timepoint)) %>%
  drop_na(ftc_mean) %>%
  ggplot(aes(x = timepoint,
             y = ftc_mean)) +
  geom_violin(aes(fill = timepoint), 
              lwd = 0.25) +
  geom_boxplot(width=0.2, lwd = 0.25, outlier.shape = NA) +
  scale_y_log10() +
  scale_fill_manual(values = c("white", "white", "white","white")) +
  theme_bw() +
  theme(axis.text.x = element_text(angle = 90, vjust = 0.5, hjust=1))+
  labs(x = "timepoint", y = 'ftc_mean')+
  theme(text = element_text(size=8), 
        axis.text	=element_text(color='black', size=8), 
        plot.title = element_text(size=8),
        legend.position = 'none') +
  facet_wrap(~ condition)

pdf("S1C_raw_conversionsGFP_dnCaf1.pdf", width = 4, height = 2)
gfp_caf
dev.off()


## load conversion rates for subcellular compartments

conversions_comp <- read.csv('../data/ftc_raw_compartments.csv') %>% 
  filter(group %in% c("Cyto", "Pro", "Nuc")) %>% 
  separate(gene, into = c("gene_id", "gene_name"), sep = "_") %>% 
  select(group, gene_id, timepoint, b.rep, cr_raw)

hl_estimated_genes <- read.csv('../data/HL.csv') %>% 
  filter(group %in% c("Cyto", "Pro", "Nuc")) %>% 
  select(gene_id, group) 

merge <- merge(conversions_comp, hl_estimated_genes, by = c('gene_id', 'group')) %>% 
  group_by(group, gene_id, timepoint) %>% 
  mutate(cr_raw_mean = mean(cr_raw, na.rm = T)) %>% 
  ungroup() %>% 
  select(-b.rep,
         -cr_raw) %>% 
  distinct()

comps <- 
  merge %>%
  mutate(timepoint = as.factor(timepoint)) %>%
  drop_na(cr_raw_mean) %>%
  ggplot(aes(x = timepoint,
             y = cr_raw_mean/100)) +
  geom_violin(aes(fill = timepoint), 
              lwd = 0.25) +
  geom_boxplot(width=0.2, lwd = 0.25, outlier.shape = NA) +
  scale_y_log10() +
  scale_fill_manual(values = c("white", "white", "white","white")) +
  theme_bw() +
  theme(axis.text.x = element_text(angle = 90, vjust = 0.5, hjust=1))+
  labs(x = "timepoint", y = 'ftc_mean')+
  theme(text = element_text(size=8), 
        axis.text	=element_text(color='black', size=8), 
        plot.title = element_text(size=8),
        legend.position = 'none') +
  facet_wrap(~ group)


# saving as pdf
pdf("S1C_raw_conversions_compartments.pdf", width = 7, height = 2)
comps
dev.off()
