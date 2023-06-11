library (tidyverse)
library (dplyr)
library (ggplot2)
library(gridExtra)
library(ggrepel)
library(rstatix)
library(ggpubr)


## load and merge HL and Deseq data

hl <- read_tsv('../data/hl_table_full_niko.tsv') %>% 
  filter(filtered == 'pseudoR2', condition == 'GFP') %>% 
  rename(hl = hl_ftc3_tc) %>% 
  select(gene_id, hl) 

deseq <- read.csv('../data/SupplementaryTable3_DE.csv') %>% 
  filter(experiment == 'PCN.ThreeWay_dnCaf1.CompSep' | experiment == 'PCN.Twoway_dnCaf1.CompSep') %>% 
  filter(comparison == 'Pro_GFP.vs.Cyto_GFP' |
           comparison == 'Neurite_GFP.vs.Soma_GFP') %>% 
  select(gene_id, log2FoldChange, comparison, DE) %>% 
  pivot_wider(values_from = log2FoldChange,
              names_from = comparison) 

merge <- merge(deseq, hl, by = 'gene_id') 

# Pearson correlation test 3C
  
Pro_GFP.vs.Cyto_GFP_list <- merge %>% 
  filter(!is.na(Pro_GFP.vs.Cyto_GFP), !is.na(hl)) %>% 
  arrange(Pro_GFP.vs.Cyto_GFP) %>% 
  select(Pro_GFP.vs.Cyto_GFP) %>% 
  unlist()

  hl_list <- merge %>% 
    filter(!is.na(Pro_GFP.vs.Cyto_GFP), !is.na(hl)) %>% 
    arrange(Pro_GFP.vs.Cyto_GFP) %>% 
    select(hl) %>% 
    unlist()
  
cor.test(Pro_GFP.vs.Cyto_GFP_list, hl_list, method = 'pearson')
  
  # p-value < 2.2e-16, cor 0.3041048 
  
  
  # Pearson correlation test 2C
  
  Neurite_GFP.vs.Soma_GFP_list <- merge %>% 
    filter(!is.na(Neurite_GFP.vs.Soma_GFP), !is.na(hl)) %>% 
    arrange(Neurite_GFP.vs.Soma_GFP) %>% 
    select(Neurite_GFP.vs.Soma_GFP) %>% 
    unlist()
  
  hl_list <- merge %>% 
    filter(!is.na(Neurite_GFP.vs.Soma_GFP), !is.na(hl)) %>% 
    arrange(Neurite_GFP.vs.Soma_GFP) %>% 
    select(hl) %>% 
    unlist()
  
cor.test(Neurite_GFP.vs.Soma_GFP_list, hl_list, method = 'pearson')
  
  # p-value < 2.2e-16, cor 0.3972613 


# plotting percentiles of HL against localization 3C

loc.hl.3c <- merge %>% 
  filter(!is.na(Pro_GFP.vs.Cyto_GFP), !is.na(hl)) %>% 
  mutate(percentile_hl = cut(hl, quantile(hl, seq(0,1, .1), na.rm = T), include_lowest = TRUE)) %>% 
  arrange(percentile_hl) %>% 
  filter(!is.na(percentile_hl)) %>% 
  ggplot(aes(percentile_hl,
                              Pro_GFP.vs.Cyto_GFP)) +
  geom_boxplot(aes(fill = percentile_hl), lwd=0.25, outlier.size = 0.1, outlier.color = 'darkgrey')+
  geom_hline(yintercept = 0) +
  theme_bw()+
  theme(axis.text.x = element_text(angle = 90, vjust = 0.5, hjust=1))+
  labs(x = "percentile HL", y = 'RNA neurite/soma (log2FC)')+
  ggtitle('Pearson 0.3041048, p < 2.2e-16')+
  scale_fill_manual(values = alpha(rep("red", 10), seq(0.2, 0.9, length.out = 10))) +
  ylim(c(-3,6)) +
  theme(text = element_text(size=8), 
        axis.text	=element_text(color='black', size=8), 
        plot.title = element_text(size=8)) +
  theme(legend.position = 'none')


loc.hl.2c <- merge %>% 
  filter(!is.na(Neurite_GFP.vs.Soma_GFP), !is.na(hl)) %>% 
  mutate(percentile_hl = cut(hl, quantile(hl, seq(0,1, .1), na.rm = T), include_lowest = TRUE)) %>% 
  arrange(percentile_hl) %>% 
  filter(!is.na(percentile_hl)) %>% 
  ggplot(aes(percentile_hl,
             Neurite_GFP.vs.Soma_GFP)) +
  geom_boxplot(aes(fill = percentile_hl), lwd=0.25, outlier.size = 0.1, outlier.color = 'darkgrey')+
  geom_hline(yintercept = 0) +
  theme_bw()+
  theme(axis.text.x = element_text(angle = 90, vjust = 0.5, hjust=1))+
  labs(x = "percentile HL", y = 'RNA neurite/soma (log2FC)')+
  ggtitle('Pearson 0.3972613, p < 2.2e-16')+
  scale_fill_manual(values = alpha(rep("red", 10), seq(0.2, 0.9, length.out = 10))) +
  ylim(c(-3,6)) +
  theme(text = element_text(size=8), 
        axis.text	=element_text(color='black', size=8), 
        plot.title = element_text(size=8)) +
  theme(legend.position = 'none')


pdf("1D_loc.vs.hl.3c.pdf", width = 3, height = 3)
loc.hl.3c
dev.off()

pdf("S1E_loc.hl.2c.pdf", width = 3, height = 3)
loc.hl.2c
dev.off()


