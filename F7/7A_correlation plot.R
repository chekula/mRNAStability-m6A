library (tidyverse)
library (dplyr)
library (ggplot2)
library(gridExtra)
library(ggrepel)
library(rstatix)
library('ggforce')


# load the data

are <- read.csv('../data/ARE.csv') %>% 
  select(gene_id, gene_name, TPU.length, ARE.cts) 

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

codonopt <- read.csv('../data/codon.opt.mouse.csv') %>% 
  select(gene_name, tAI) 

# merge the data
data <- merge(m6a.sysy, are, by = 'gene_id')
data <- merge(data, codonopt, by = 'gene_name')
data <- merge(hl, data, by = 'gene_id') %>% 
  drop_na()

#--------- correlation plot

# Install and load the corrplot package
library(corrplot)

# Create a matrix of correlation coefficients
corplot.data <- data[, c(2, 4, 6, 7)] 
cor_matrix <- cor(corplot.data)

write.csv(cor_matrix, '7A_cor_matrix.csv')

# Plot the correlation matrix using the corrplot function

pdf("7A_corrplot_colors_reversed.pdf", width = 5, height = 5)
corrplot(cor_matrix, method = 'ellipse', col = rev(colorRampPalette(c("#67001F", "#B2182B", "#D6604D", "#F4A582", "#FDDBC7", "#FFFFFF", "#D1E5F0", "#92C5DE", "#4393C3", "#2166AC", "#053061"))(10)))
dev.off()
