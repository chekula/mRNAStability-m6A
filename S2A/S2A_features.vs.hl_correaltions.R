library (tidyverse)
library (dplyr)
library (ggplot2)
library(gridExtra)
library(ggrepel)
library(rstatix)
library('ggforce')


# load the data

hl <- read_tsv('../data/hl_table_full_niko.tsv') %>% 
  filter(filtered == 'pseudoR2') %>% 
  select(gene_id,
         condition,
         hl_ftc3_tc) %>% 
  rename(hl = hl_ftc3_tc) %>% 
  pivot_wider(names_from = condition,
              values_from = hl,
              names_prefix = 'hl.')

are <- read.csv('../data/ARE.csv') %>% 
  select(gene_id, gene_name, TPU.length, ARE.cts) 

m6a.sysy <- read.csv('../data/m6A.csv') %>% 
  filter(dataset_msa == 'sysy_total') %>% 
  select(gene_id, log2FoldChange) %>% 
  rename(m6a.enrich = log2FoldChange) 

codonopt <- read.csv('../data/codon.opt.mouse.csv') %>% 
  select(gene_name, tAI) 

data <- merge(m6a.sysy, are, by = 'gene_id')
data <- merge(data, codonopt, by = 'gene_name')
data <- merge(hl, data, by = 'gene_id', all = TRUE)

write.csv(data, '../data/features_combined.csv', row.names = F)


#--------- correlation plot

# Install and load the corrplot package
library(corrplot)
library(Hmisc)

# Select the columns of interest and drop missing values
corplot.data <- data[, c(2, 3, 5, 6, 7, 8)] %>% drop_na()

# Calculate the correlation matrix and associated p-values
cor_matrix <- rcorr(as.matrix(corplot.data), type = "pearson")

# Extract the correlation coefficients and associated p-values
R <- cor_matrix$r
p <- cor_matrix$P

# Flatten the p-value matrix into a vector
p_vector <- as.vector(p)

# Adjust p-values for multiple testing using  Benjamini-Hochberg procedure
p.adjusted_vector <- p.adjust(p_vector, method = "BH")

# Reshape the adjusted p-values back into a matrix of the original shape
p.adjusted <- matrix(p.adjusted_vector, nrow = nrow(p), ncol = ncol(p))

# Preserve the names of the rows and columns
rownames(p.adjusted) <- rownames(p)
colnames(p.adjusted) <- colnames(p)

write.csv(R, 'correlation_R_S2A.csv')
write.csv(p, 'correlation_p_S2A.csv')
write.csv(p.adjusted, 'correlation_p_adjusted_S2A.csv')

# Plotting S2A: m6A

m6a.gfp <- 
  data[, c(2, 3, 5, 6, 7, 8)] %>% 
  drop_na() %>% 
  mutate(percentile_m6a = cut(m6a.enrich,
                              quantile(m6a.enrich, seq(0,1, .25), 
                                       na.rm = T), 
                              include_lowest = TRUE)) %>% 
  filter(!is.na(percentile_m6a)) %>% 
  ggplot(aes(x = percentile_m6a,
             y = hl.GFP,
             fill = percentile_m6a)) +
  geom_boxplot(lwd = 0.15, outlier.shape = NA) +
  scale_fill_manual(values = alpha(rep("red", 4), seq(0.2, 0.8, length.out = 4))) +
  theme_bw() +
  theme(axis.text.x = element_text(angle = 90, vjust = 0.5, hjust=1),
        axis.text.y = element_text(size=5),
        axis.ticks = element_line(size = 0.1),
        axis.line = element_line(size = 0.1),
        plot.title = element_text(size = 5),
        legend.position = 'none') +
  ylim(c(0, 14)) +
  ggtitle('Pearson R = -0.4182756,
          p < 2.2e-16')

m6a.caf <- 
  data[, c(2, 3, 5, 6, 7, 8)] %>% 
  drop_na() %>% 
  mutate(percentile_m6a = cut(m6a.enrich,
                              quantile(m6a.enrich, seq(0,1, .25), 
                                       na.rm = T), 
                              include_lowest = TRUE)) %>% 
  filter(!is.na(percentile_m6a)) %>% 
  ggplot(aes(x = percentile_m6a,
             y = hl.dnCaf1,
             fill = percentile_m6a)) +
  geom_boxplot(lwd = 0.15, outlier.shape = NA) +
  scale_fill_manual(values = alpha(rep("red", 4), seq(0.2, 0.8, length.out = 4))) +
  theme_bw() +
  theme(axis.text.x = element_text(angle = 90, vjust = 0.5, hjust=1),
        axis.text.y = element_text(size=5),
        axis.ticks = element_line(size = 0.1),
        axis.line = element_line(size = 0.1),
        plot.title = element_text(size = 5),
        legend.position = 'none') +
  ylim(c(0, 14)) +
  ggtitle('Pearson R = -0.3832778,
          p < 2.2e-16')

pdf("S2A_m6a.pdf", width = 3, height = 3)
grid.arrange(m6a.gfp, m6a.caf, ncol= 2)
dev.off()

# Plotting S2A: ARE

are.gfp <- 
  data[, c(2, 3, 5, 6, 7, 8)] %>% 
  drop_na() %>% 
  mutate(percentile_ARE = cut(ARE.cts,
                              quantile(ARE.cts, seq(0,1, .25), 
                                       na.rm = T), 
                              include_lowest = TRUE)) %>% 
  filter(!is.na(percentile_ARE)) %>% 
  ggplot(aes(x = percentile_ARE,
             y = hl.GFP,
             fill = percentile_ARE)) +
  geom_boxplot(lwd = 0.15, outlier.shape = NA) +
  scale_fill_manual(values = alpha(rep("red", 4), seq(0.2, 0.8, length.out = 4))) +
  theme_bw() +
  theme(axis.text.x = element_text(angle = 90, vjust = 0.5, hjust=1),
        axis.text.y = element_text(size=5),
        axis.ticks = element_line(size = 0.1),
        axis.line = element_line(size = 0.1),
        plot.title = element_text(size = 5),
        legend.position = 'none') +
  ylim(c(0, 14)) +
  ggtitle('Pearson R = -0.1648989,
          p < 2.2e-16')

are.caf <- 
  data[, c(2, 3, 5, 6, 7, 8)] %>% 
  drop_na() %>% 
  mutate(percentile_ARE = cut(ARE.cts,
                              quantile(ARE.cts, seq(0,1, .25), 
                                       na.rm = T), 
                              include_lowest = TRUE)) %>% 
  filter(!is.na(percentile_ARE)) %>% 
  ggplot(aes(x = percentile_ARE,
             y = hl.dnCaf1,
             fill = percentile_ARE)) +
  geom_boxplot(lwd = 0.15, outlier.shape = NA) +
  scale_fill_manual(values = alpha(rep("red", 4), seq(0.2, 0.8, length.out = 4))) +
  theme_bw() +
  theme(axis.text.x = element_text(angle = 90, vjust = 0.5, hjust=1),
        axis.text.y = element_text(size=5),
        axis.ticks = element_line(size = 0.1),
        axis.line = element_line(size = 0.1),
        plot.title = element_text(size = 5),
        legend.position = 'none') +
  ylim(c(0, 14)) +
  ggtitle('Pearson R = -0.1712564,
          p < 2.2e-16')

pdf("S2A_are.pdf", width = 3, height = 3)
grid.arrange(are.gfp, are.caf, ncol= 2)
dev.off()

# Plotting S2A: 3UTR

tpu.gfp <- 
  data[, c(2, 3, 5, 6, 7, 8)] %>% 
  drop_na() %>% 
  mutate(percentile_TPU = cut(TPU.length,
                              quantile(TPU.length, seq(0,1, .25), 
                                       na.rm = T), 
                              include_lowest = TRUE)) %>% 
  filter(!is.na(percentile_TPU)) %>% 
  ggplot(aes(x = percentile_TPU,
             y = hl.GFP,
             fill = percentile_TPU)) +
  geom_boxplot(lwd = 0.15, outlier.shape = NA) +
  scale_fill_manual(values = alpha(rep("red", 4), seq(0.2, 0.8, length.out = 4))) +
  theme_bw() +
  theme(axis.text.x = element_text(angle = 90, vjust = 0.5, hjust=1),
        axis.text.y = element_text(size=5),
        axis.ticks = element_line(size = 0.1),
        axis.line = element_line(size = 0.1),
        plot.title = element_text(size = 5),
        legend.position = 'none') +
  ylim(c(0, 14)) +
  ggtitle('Pearson R = -0.1298628,
          p < 2.2e-16')

tpu.caf <- 
  data[, c(2, 3, 5, 6, 7, 8)] %>% 
  drop_na() %>% 
  mutate(percentile_TPU = cut(TPU.length,
                              quantile(TPU.length, seq(0,1, .25), 
                                       na.rm = T), 
                              include_lowest = TRUE)) %>% 
  filter(!is.na(percentile_TPU)) %>% 
  ggplot(aes(x = percentile_TPU,
             y = hl.dnCaf1,
             fill = percentile_TPU)) +
  geom_boxplot(lwd = 0.15, outlier.shape = NA) +
  scale_fill_manual(values = alpha(rep("red", 4), seq(0.2, 0.8, length.out = 4))) +
  theme_bw() +
  theme(axis.text.x = element_text(angle = 90, vjust = 0.5, hjust=1),
        axis.text.y = element_text(size=5),
        axis.ticks = element_line(size = 0.1),
        axis.line = element_line(size = 0.1),
        plot.title = element_text(size = 5),
        legend.position = 'none') +
  ylim(c(0, 14)) +
  ggtitle('Pearson R = -0.1514933,
          p < 2.2e-16')

pdf("S2A_tpu.pdf", width = 3, height = 3)
grid.arrange(tpu.gfp, tpu.caf, ncol= 2)
dev.off()

# Plotting S2A: codon optimality

tai.gfp <- 
  data[, c(2, 3, 5, 6, 7, 8)] %>% 
  drop_na() %>% 
  mutate(percentile_tai = cut(tAI,
                              quantile(tAI, seq(0,1, .25), 
                                       na.rm = T), 
                              include_lowest = TRUE)) %>% 
  filter(!is.na(percentile_tai)) %>% 
  ggplot(aes(x = percentile_tai,
             y = hl.GFP,
             fill = percentile_tai)) +
  geom_boxplot(lwd = 0.15, outlier.shape = NA) +
  scale_fill_manual(values = alpha(rep("red", 4), seq(0.2, 0.8, length.out = 4))) +
  theme_bw() +
  theme(axis.text.x = element_text(angle = 90, vjust = 0.5, hjust=1),
        axis.text.y = element_text(size=5),
        axis.ticks = element_line(size = 0.1),
        axis.line = element_line(size = 0.1),
        plot.title = element_text(size = 5),
        legend.position = 'none') +
  ylim(c(0, 14)) +
  ggtitle('Pearson R = 0.1697477,
          p < 2.2e-16')

tai.caf <- 
  data[, c(2, 3, 5, 6, 7, 8)] %>% 
  drop_na() %>% 
  mutate(percentile_tai = cut(tAI,
                              quantile(tAI, seq(0,1, .25), 
                                       na.rm = T), 
                              include_lowest = TRUE)) %>% 
  filter(!is.na(percentile_tai)) %>% 
  ggplot(aes(x = percentile_tai,
             y = hl.dnCaf1,
             fill = percentile_tai)) +
  geom_boxplot(lwd = 0.15, outlier.shape = NA) +
  scale_fill_manual(values = alpha(rep("red", 4), seq(0.2, 0.8, length.out = 4))) +
  theme_bw() +
  theme(axis.text.x = element_text(angle = 90, vjust = 0.5, hjust=1),
        axis.text.y = element_text(size=5),
        axis.ticks = element_line(size = 0.1),
        axis.line = element_line(size = 0.1),
        plot.title = element_text(size = 5),
        legend.position = 'none') +
  ylim(c(0, 14)) +
  ggtitle('Pearson R = 0.1634230,
          p < 2.2e-16')

pdf("S2A_tAI.pdf", width = 3, height = 3)
grid.arrange(tai.gfp, tai.caf, ncol= 2)
dev.off()

