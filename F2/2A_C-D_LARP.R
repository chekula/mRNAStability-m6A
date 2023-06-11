library (tidyverse)
library (dplyr)
library (ggplot2)
library(gridExtra)
library(rstatix)
library(ggpubr)
library("cowplot")
library('ggforce')


## load data

deseq.scr <- read_tsv('../data/scr_neu.vs.scr_cyto.deseq_results.tsv') %>% 
  select(gene_id, log2FoldChange, geneName) %>% 
  rename(scr_neu.vs.scr_cyto = log2FoldChange, gene_name = geneName) 

deseq.larp <- read_tsv('../data/larp_neu.vs.larp_cyto.deseq_results.tsv') %>% 
  select(gene_id, log2FoldChange) %>% 
  rename(larp_neu.vs.larp_cyto = log2FoldChange) 

deseq.neu <- read_tsv('../data/larp_neu.vs.scr_neu.deseq_results.tsv') %>% 
  select(gene_id, log2FoldChange, baseMean) %>% 
  rename(larp_neu.vs.scr_neu = log2FoldChange,
         baseMean_neu = baseMean) 

deseq.cyto <- read_tsv('../data/larp_cyto.vs.scr_cyto.deseq_results.tsv') %>% 
  select(gene_id, log2FoldChange, baseMean) %>% 
  rename(larp_cyto.vs.scr_cyto = log2FoldChange,
         baseMean_cyto = baseMean) 


data <- merge(deseq.scr, deseq.larp, by = 'gene_id') 
data <- merge(data, deseq.neu, by = 'gene_id') 
data <- merge(data, deseq.cyto, by = 'gene_id') 
data <- data %>% 
  mutate(delta_larp.vs.scr = larp_neu.vs.larp_cyto - scr_neu.vs.scr_cyto)

# load 5TOP gene names
top_genes <- read.csv("../data/TOP_gene_names.csv")

# Filter the rows with gene names starting with "Rp"
rp_genes <- top_genes %>% 
  filter(str_starts(gene_name, "Rp")) %>% 
  pull(gene_name)

# Mutate the data and add the group column
data5top <- data %>% 
  mutate(group = if_else(gene_name %in% rp_genes, "RP.5TOP", "else"))

# subsetting a table with 'RP.5TOP' only
toponly <- data5top %>% 
  filter(group == 'RP.5TOP')

## - Welch variation of t-test: is for unequal variance in compared samples: localization between RP.TOP and else

rp5top <- data5top %>% 
  filter(group == 'RP.5TOP')

rest <- data5top %>% 
  filter(group == 'else')

statistic1 <- 
  t.test(rp5top$scr_neu.vs.scr_cyto, rest$scr_neu.vs.scr_cyto, paired=F, var.equal = F) 


# plot changes in localization Larp KD vs scrambled

p1 <- 
  ggplot(data5top, aes(x = group,
             y = scr_neu.vs.scr_cyto)) +
  geom_boxplot(aes(fill = group), alpha = 0.5, lwd = 0.15, outlier.shape = NA)+
  theme(text = element_text(size=12), axis.text	=element_text(color='black', size=12), plot.title = element_text(size=2)) +
  geom_hline(yintercept = 0, lwd=0.15) +
  theme_bw()+
  theme(axis.text.x = element_text(angle = 90, vjust = 0.5, hjust=1))+
  labs(x = "localization neu.vs.cyto")+
  scale_color_manual(values = c('darkgrey', 'red'))+
  scale_fill_manual(values = c('darkgrey', 'red'))+
  ylim(c(-4,5)) +
  theme(axis.text.x = element_text(angle = 90, vjust = 0.5, hjust=1),
        axis.text.y = element_text(size=5),
        axis.ticks = element_line(size = 0.1),
        axis.line = element_line(size = 0.1),
        plot.title = element_text(size = 5),
        legend.position = 'none') +
  geom_sina(toponly, mapping = aes(x = group,
                                   y = scr_neu.vs.scr_cyto),
            color = 'red',
            size = 0.1,
            alpha = 0.5) +
  ggtitle(paste('Welch p-value \n', statistic1[3]))


statistic2 <- 
  t.test(rp5top$delta_larp.vs.scr, rest$delta_larp.vs.scr, paired=F, var.equal = F) 

p2 <- 
    ggplot(data5top, aes(x = group,
             y = delta_larp.vs.scr)) +
  geom_boxplot(aes(fill = group), alpha = 0.5, lwd = 0.15, outlier.shape = NA)+
  theme(text = element_text(size=12), axis.text	=element_text(color='black', size=12), plot.title = element_text(size=2)) +
  geom_hline(yintercept = 0, lwd=0.15) +
  theme_bw()+
  theme(axis.text.x = element_text(angle = 90, vjust = 0.5, hjust=1))+
  ylim(c(-4,5)) +
  labs(x = "delta localization larp-scr") +
  scale_color_manual(values = c('darkgrey', 'red'))+
  scale_fill_manual(values = c('darkgrey', 'red'))+
  theme(axis.text.x = element_text(angle = 90, vjust = 0.5, hjust=1),
        axis.text.y = element_text(size=5),
        axis.ticks = element_line(size = 0.1),
        axis.line = element_line(size = 0.1),
        plot.title = element_text(size = 5),
        legend.position = 'none') +
    geom_sina(toponly, mapping = aes(x = group,
                                     y = delta_larp.vs.scr),
              color = 'red',
              size = .1,
              alpha = 0.5) +
  ggtitle(paste('Welch p-value \n', statistic2[3]))

## saving as pdf
pdf("2A_D_5top.loc.pdf", width = 3, height = 3)
grid.arrange(p1, p2, ncol = 2)
dev.off()


# boxplot changes in RNA levels larp KD vs scrambled

statistic3 <- 
  t.test(rp5top$larp_neu.vs.scr_neu, rest$larp_neu.vs.scr_neu, paired=F, var.equal = F) 

p3 <- 
    ggplot(data5top, aes(x = group,
                         y = larp_neu.vs.scr_neu)) +
    geom_boxplot(aes(fill = group), alpha = 0.5, lwd=0.15, outlier.shape = NA)+
    theme(text = element_text(size=12), axis.text	=element_text(color='black', size=12), plot.title = element_text(size=2)) +
    geom_hline(yintercept = 0, lwd=0.15) +
    theme_bw()+
    theme(axis.text.x = element_text(angle = 90, vjust = 0.5, hjust=1)) +
  scale_fill_manual(values = c('darkgrey', 'red'))+
  theme(axis.text.x = element_text(angle = 90, vjust = 0.5, hjust=1),
        axis.text.y = element_text(size=5),
        axis.ticks = element_line(size = 0.1),
        axis.line = element_line(size = 0.1),
        plot.title = element_text(size = 5),
        legend.position = 'none') +
geom_sina(toponly, mapping = aes(x = group,
                                 y = larp_neu.vs.scr_neu),
          color = 'red',
          size = .1,
          alpha = 0.5) +
  ggtitle(paste('Welch p-value \n', statistic3[3]))+
  ylim(c(-4,3))


statistic4 <- 
  t.test(rp5top$larp_cyto.vs.scr_cyto, rest$larp_cyto.vs.scr_cyto, paired=F, var.equal = F) 


p4 <- 
  ggplot(data5top, aes(x = group,
                       y = larp_cyto.vs.scr_cyto)) +
  geom_boxplot(aes(fill = group), alpha = 0.5, lwd=0.15, outlier.shape = NA)+
  theme(text = element_text(size=12), axis.text	=element_text(color='black', size=12), plot.title = element_text(size=2)) +
  geom_hline(yintercept = 0, lwd=0.15) +
  theme_bw()+
  scale_fill_manual(values = c('darkgrey', 'red'))+
  theme(axis.text.x = element_text(angle = 90, vjust = 0.5, hjust=1),
        axis.text.y = element_text(size=5),
        axis.ticks = element_line(size = 0.1),
        axis.line = element_line(size = 0.1),
        plot.title = element_text(size = 5),
        legend.position = 'none') +
  geom_sina(toponly, mapping = aes(x = group,
                                   y = larp_cyto.vs.scr_cyto),
            color = 'red',
            size = .1,
            alpha = 0.5) +
  ggtitle(paste('Welch p-value \n', statistic4[3]))+
  ylim(c(-4,3))

## saving as pdf
pdf("2C_5top.RNAlevels.pdf", width = 3, height = 3)
grid.arrange(p3, p4, ncol = 2)
dev.off()




