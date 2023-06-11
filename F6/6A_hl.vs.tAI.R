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
         gene_id,
         gene_name) %>% 
  rename(hl = hl_ftc3_tc)

codonopt <- read.csv('../data/codon.opt.mouse.csv') %>% 
  select(gene_name, tAI) 

data <- merge(codonopt, hl, by = 'gene_name')

# add old hl for mESc-derived soma and neurotes

old.hl <- read.csv('../data/HL.csv') %>%
  select(gene_id, group, Half_life) %>% 
  filter(group %in% c('Neurite', 'Soma')) %>%
  pivot_wider(names_from = group, values_from = Half_life, names_prefix = "hl.")

combined.data <- merge(data, old.hl, by = 'gene_id', all = T)

#---------------------- hl vs. codon optimality

# Pearson correlation test PCN

hl_list <- combined.data %>% 
  filter(!is.na(hl), !is.na(tAI)) %>% 
  arrange(hl)  %>% 
  select(hl) %>% 
  unlist()

tAI_list <- combined.data %>% 
  filter(!is.na(hl), !is.na(tAI)) %>% 
  arrange(hl) %>% 
  select(tAI) %>% 
  unlist()

Pearson <- cor.test(tAI_list, hl_list, method = 'pearson') 

# Pearson cor 0.1780812, p-value < 2.2e-16

# plotting
hl.tAI.pcn <- 
  combined.data %>%
  filter(!is.na(hl)) %>%
  mutate(percentile = cut(hl, quantile(hl, seq(0,1, .2), na.rm = T), include_lowest = TRUE)) %>%
  filter(!is.na(percentile)) %>%
  arrange(percentile) %>%
  ggplot(aes(percentile, tAI, fill = percentile)) +
  geom_boxplot(lwd = 0.15, outlier.shape = NA) +
  scale_fill_manual(values = alpha(rep("red", 5), seq(0.2, 0.9, length.out = 5))) +
  theme_bw() +
  theme(axis.text.x = element_text(angle = 90, vjust = 0.5, hjust=1),
        axis.text.y = element_text(size=5),
        axis.ticks = element_line(size = 0.1),
        axis.line = element_line(size = 0.1),
        plot.title = element_text(size = 5),
        legend.position = "none") +
  ggtitle('PCN, Pearson cor 0.1780812 \n p < 2.2e-16') +
  labs(x = 'hl percentile', y = 'tAI') +
  ylim(c(0.17,0.23))


# Pearson correlation test Soma mESC

hl_list <- combined.data %>% 
  filter(!is.na(hl.Soma), !is.na(tAI)) %>% 
  arrange(hl.Soma) %>% 
  select(hl.Soma) %>% 
  unlist()

tAI_list <- combined.data %>% 
  filter(!is.na(hl.Soma), !is.na(tAI)) %>% 
  arrange(hl.Soma) %>% 
  select(tAI) %>% 
  unlist()

Pearson <- cor.test(tAI_list, hl_list, method = 'pearson') 

# Pearson cor 0.09522988  p-value = 4.834e-05

# plotting
hl.tAI.soma.mESC <- 
  combined.data %>%
  filter(!is.na(hl.Soma)) %>%
  mutate(percentile = cut(hl.Soma, quantile(hl.Soma, seq(0,1, .2), na.rm = T), include_lowest = TRUE)) %>%
  filter(!is.na(percentile)) %>%
  arrange(percentile) %>%
  ggplot(aes(percentile, tAI, fill = percentile)) +
  geom_boxplot(lwd = 0.15, outlier.shape = NA) +
  scale_fill_manual(values = alpha(rep("red", 5), seq(0.2, 0.9, length.out = 5))) +
  theme_bw() +
  theme(axis.text.x = element_text(angle = 90, vjust = 0.5, hjust=1),
        axis.text.y = element_text(size=5),
        axis.ticks = element_line(size = 0.1),
        axis.line = element_line(size = 0.1),
        plot.title = element_text(size = 5),
        legend.position = "none") +
  ggtitle('PCN,Pearson cor 0.09522988 \n p = 4.834e-05') +
  ylim(c(0.17, 0.23))+
  labs(x = 'hl percentile, mESC Soma', y = 'tAI')

pdf("6A_hl.tAI.pdf", width = 3, height = 3)
grid.arrange(hl.tAI.pcn, hl.tAI.soma.mESC, ncol = 2)
dev.off()


# Pearson correlation test neurites mESC

hl_list <- combined.data %>% 
  filter(!is.na(hl.Neurite), !is.na(tAI)) %>% 
  arrange(hl.Neurite) %>% 
  select(hl.Neurite) %>% 
  unlist()

tAI_list <- combined.data %>% 
  filter(!is.na(hl.Neurite), !is.na(tAI)) %>% 
  arrange(hl.Neurite) %>% 
  select(tAI) %>% 
  unlist()

Pearson <- cor.test(tAI_list, hl_list, method = 'pearson') 

# Pearson cor 0.07295988  p-value = 0.001869

# plotting
hl.tAI.neurite.mESC <- 
  combined.data %>%
  filter(!is.na(hl.Neurite)) %>%
  mutate(percentile = cut(hl.Neurite, quantile(hl.Neurite, seq(0,1, .2), na.rm = T), include_lowest = TRUE)) %>%
  filter(!is.na(percentile)) %>%
  arrange(percentile) %>%
  ggplot(aes(percentile, tAI, fill = percentile)) +
  geom_boxplot(lwd = 0.15, outlier.shape = NA) +
  scale_fill_manual(values = alpha(rep("red", 5), seq(0.2, 0.9, length.out = 5))) +
  theme_bw() +
  theme(axis.text.x = element_text(angle = 90, vjust = 0.5, hjust=1),
        axis.text.y = element_text(size=5),
        axis.ticks = element_line(size = 0.1),
        axis.line = element_line(size = 0.1),
        plot.title = element_text(size = 5),
        legend.position = "none") +
  ggtitle('PCN,Pearson cor 0.07295988 \n p = 0.001869') +
  ylim(c(0.17, 0.23))+
  labs(x = 'hl percentile, mESC Neurites', y = 'tAI')

