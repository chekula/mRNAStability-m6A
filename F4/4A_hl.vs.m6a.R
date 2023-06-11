library (tidyverse)
library (dplyr)
library (ggplot2)
library(gridExtra)
library(ggrepel)
library(rstatix)
library('ggforce')


# load the data

hl <- read_tsv('../data/hl_table_full_niko.tsv') %>% 
  filter(condition == "GFP", 
         filtered == 'pseudoR2') %>% 
  select(hl_ftc3_tc,
         gene_id) %>% 
  rename(hl = hl_ftc3_tc)

m6a.sysy <- read.csv('../data/m6A.csv') %>% 
  filter(dataset_msa == 'sysy_total') %>% 
  select(gene_id, log2FoldChange) %>% 
  rename(m6a.enrich = log2FoldChange) 

data <- merge(m6a.sysy, hl, by = 'gene_id')

# add old hl for mESc-derived soma nad neurotes

old.hl <- read.csv('../data/HL.csv') %>%
  select(gene_id, group, Half_life) %>% 
  filter(group %in% c('Neurite', 'Soma')) %>%
  pivot_wider(names_from = group, values_from = Half_life, names_prefix = "hl.")

combined.data <- merge(data, old.hl, by = 'gene_id', all = T)


# Pearson correlation test PCN

hl_list <- combined.data %>% 
  filter(!is.na(hl), !is.na(m6a.enrich)) %>% 
  arrange(hl) %>% 
  select(hl) %>% 
  unlist()

m6a_list <- combined.data %>% 
  filter(!is.na(hl), !is.na(m6a.enrich)) %>% 
  arrange(hl) %>% 
  select(m6a.enrich) %>% 
  unlist()

Pearson <- cor.test(m6a_list, hl_list, method = 'pearson') 

# Pearson cor -0.436742, p-value < 2.2e-16

# plotting
hl.m6a.pcn <- 
  combined.data %>%
  filter(!is.na(hl)) %>%
  mutate(percentile = cut(hl, quantile(hl, seq(0,1, .1), na.rm = T), include_lowest = TRUE)) %>%
  filter(!is.na(percentile)) %>%
  arrange(percentile) %>%
  ggplot(aes(percentile, m6a.enrich, fill = percentile)) +
  geom_boxplot(lwd = 0.15, outlier.shape = NA) +
  scale_fill_manual(values = alpha(rep("red", 10), seq(0.2, 0.9, length.out = 10))) +
  theme_bw() +
  theme(axis.text.x = element_text(angle = 90, vjust = 0.5, hjust=1),
        axis.text.y = element_text(size=5),
        axis.ticks = element_line(size = 0.1),
        axis.line = element_line(size = 0.1),
        plot.title = element_text(size = 5),
        legend.position = "none") +
  ggtitle('PCN,Pearson cor -0.436742 \n p < 2.2e-16') +
  ylim(c(-6, 6))+
  labs(x = 'hl percentile', y = 'm6a enrichment sysy')



# Pearson correlation test Soma mESC

hl_list <- combined.data %>% 
  filter(!is.na(hl.Soma), !is.na(m6a.enrich)) %>% 
  arrange(hl.Soma) %>% 
  select(hl.Soma) %>% 
  unlist()

m6a_list <- combined.data %>% 
  filter(!is.na(hl.Soma), !is.na(m6a.enrich)) %>% 
  arrange(hl.Soma) %>% 
  select(m6a.enrich) %>% 
  unlist()

Pearson <- cor.test(m6a_list, hl_list, method = 'pearson') 

# Pearson cor -0.3319755 p-value < 2.2e-16

# plotting
hl.m6a.soma.mESC <- 
  combined.data %>%
  filter(!is.na(hl.Soma)) %>%
  mutate(percentile = cut(hl.Soma, quantile(hl.Soma, seq(0,1, .1), na.rm = T), include_lowest = TRUE)) %>%
  filter(!is.na(percentile)) %>%
  arrange(percentile) %>%
  ggplot(aes(percentile, m6a.enrich, fill = percentile)) +
  geom_boxplot(lwd = 0.15, outlier.shape = NA) +
  scale_fill_manual(values = alpha(rep("red", 10), seq(0.2, 0.9, length.out = 10))) +
  theme_bw() +
  theme(axis.text.x = element_text(angle = 90, vjust = 0.5, hjust=1),
        axis.text.y = element_text(size=5),
        axis.ticks = element_line(size = 0.1),
        axis.line = element_line(size = 0.1),
        plot.title = element_text(size = 5),
        legend.position = "none") +
  ggtitle('mESC soma, Pearson cor -0.3319755 \n p < 2.2e-16') +
  ylim(c(-6, 6))+
  labs(x = 'hl percentile, mESC Soma', y = 'm6a enrichment sysy')

pdf("4A_hl.vs.m6a.pdf", width = 3, height = 3)
grid.arrange(hl.m6a.pcn, hl.m6a.soma.mESC, ncol = 2)
dev.off()


# same for mESC neurites (not included)

# Pearson correlation test neurites mESC

hl_list <- combined.data %>% 
  filter(!is.na(hl.Neurite), !is.na(m6a.enrich)) %>% 
  arrange(hl.Neurite) %>% 
  select(hl.Neurite) %>% 
  unlist()

m6a_list <- combined.data %>% 
  filter(!is.na(hl.Neurite), !is.na(m6a.enrich)) %>% 
  arrange(hl.Neurite)  %>% 
  select(m6a.enrich) %>% 
  unlist()

Pearson <- cor.test(m6a_list, hl_list, method = 'pearson') 

# Pearson cor -0.3099492 p-value < 2.2e-16

# plotting
hl.m6a.neurite.mESC <- 
  combined.data %>%
  filter(!is.na(hl.Neurite)) %>%
  mutate(percentile = cut(hl.Neurite, quantile(hl.Neurite, seq(0,1, .1), na.rm = T), include_lowest = TRUE)) %>%
  filter(!is.na(percentile)) %>%
  arrange(percentile) %>%
  ggplot(aes(percentile, m6a.enrich, fill = percentile)) +
  geom_boxplot(lwd = 0.15, outlier.shape = NA) +
  scale_fill_manual(values = alpha(rep("red", 10), seq(0.2, 0.9, length.out = 10))) +
  theme_bw() +
  theme(axis.text.x = element_text(angle = 90, vjust = 0.5, hjust=1),
        axis.text.y = element_text(size=5),
        axis.ticks = element_line(size = 0.1),
        axis.line = element_line(size = 0.1),
        plot.title = element_text(size = 5),
        legend.position = "none") +
  ggtitle('mESc neurite, Pearson cor -0.3099492 \n p < 2.2e-16') +
  ylim(c(-6, 6))+
  labs(x = 'hl percentile, mESC Neurites', y = 'm6a enrichment sysy')
