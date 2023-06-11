library (tidyverse)
library (dplyr)
library (ggplot2)
library(gridExtra)
library(rstatix)
library('ggforce')


# load the data

hl <- read_tsv('../data/hl_table_full_niko.tsv') %>% 
  filter(filtered == 'pseudoR2') %>% 
  select(gene_id,
         hl_ftc3_tc,
         condition) %>% 
  pivot_wider(values_from = hl_ftc3_tc,
              names_from = condition, names_prefix = "hl.")

m6a.sysy <- read.csv('../data/m6A.csv') %>% 
  filter(dataset_msa == 'sysy_total') %>% 
  select(gene_id, log2FoldChange) %>% 
  rename(m6a.enrich = log2FoldChange) 

are <- read.csv('../data/ARE.csv') %>% 
  select(gene_id, gene_name, TPU.length, ARE.cts)

codonopt <- read.csv('../data/codon.opt.mouse.csv') %>% 
  select(gene_name, tAI)

data <- merge(m6a.sysy, are, by = 'gene_id')
data <- merge(data, codonopt, by = 'gene_name')
data <- merge(hl, data, by = 'gene_id') %>% drop_na()

# add nELavl targets based on downregulation in shnElavl 

elavl.down <- read.csv('../data/SupplementaryTable3_DE.csv') %>% 
  filter(experiment == "PCN.TwoWay_shElavl234",
         comparison %in% c("Neurite_shElavl234.vs.Neurite_Scrambled",
                           "Soma_shElavl234.vs.Soma_Scrambled")) %>% 
  select(gene_id, comparison, DE) %>% 
  pivot_wider(names_from = comparison, values_from = DE) 

elavl_targets <- elavl.down %>% 
  filter(Neurite_shElavl234.vs.Neurite_Scrambled == "Scrambled" &
           Soma_shElavl234.vs.Soma_Scrambled == "Scrambled") %>% 
  pull(gene_id)

data <- data %>% 
mutate(elavl = if_else(gene_id %in% elavl_targets, "elavl.target", "none")) %>% 
 mutate(percentile_tpu = cut(TPU.length,
                             quantile(TPU.length, seq(0,1, .25), 
                                      na.rm = T), 
                             include_lowest = TRUE)) %>% 
  mutate(percentile_are = cut(ARE.cts,
                              quantile(ARE.cts, seq(0,1, .25), 
                                       na.rm = T), 
                              include_lowest = TRUE)) %>% 
  filter(!is.na(percentile_tpu)) %>% 
  filter(!is.na(percentile_are))

# ----- S6A are vs. tpu vs. hl.GFP
# Pearson correlation
unique(data$percentile_tpu)

# Define a vector to hold the p-values
p_values <- c()

# Perform correlation tests and store the p-values
low25 <- data %>% filter(percentile_tpu == '(18,698]')
hl_list <- low25 %>% select(hl.GFP) %>% unlist()
ARE_list <- low25 %>% select(ARE.cts) %>% unlist()
p_values <- c(p_values, cor.test(hl_list, ARE_list, method = 'pearson')$p.value) # -0.1232038 p-value = 5.458e-05

low50 <- data %>% filter(percentile_tpu == '(698,1.5e+03]')
hl_list <- low50 %>% select(hl.GFP) %>% unlist()
ARE_list <- low50 %>% select(ARE.cts) %>% unlist()
p_values <- c(p_values, cor.test(hl_list, ARE_list, method = 'pearson')$p.value) # -0.1395903 p-value = 7.121e-09

low75 <- data %>% filter(percentile_tpu == '(1.5e+03,2.7e+03]')
hl_list <- low75 %>% select(hl.GFP) %>% unlist()
ARE_list <- low75 %>% select(ARE.cts) %>% unlist()
p_values <- c(p_values, cor.test(hl_list, ARE_list, method = 'pearson')$p.value) # -0.1105765 p-value = 9.41e-07

low100 <- data %>% filter(percentile_tpu == '(2.7e+03,3.94e+04]')
hl_list <- low100 %>% select(hl.GFP) %>% unlist()
ARE_list <- low100 %>% select(ARE.cts) %>% unlist()
p_values <- c(p_values, cor.test(hl_list, ARE_list, method = 'pearson')$p.value) # -0.1506148 p-value = 5.602e-12

# Adjust the p-values using the Benjamini-Hochberg procedure
p_values_adjusted <- p.adjust(p_values, method = "BH")
print(p_values_adjusted)


# plot hl for ARE and TPU percentiles

are.tpu.hl <- 
  data %>% 
  ggplot(aes(x = hl.GFP,
             y = percentile_tpu,
             fill = percentile_are)) +
  geom_boxplot(lwd = 0.15, outlier.shape = NA) +
  scale_fill_manual(values = alpha(rep("red", 4), seq(0.2, 0.8, length.out = 4))) +
  theme_bw() +
  theme(axis.text.x = element_text(angle = 0, vjust = 0.5, hjust=1),
        axis.text.y = element_text(size=5),
        axis.ticks = element_line(size = 0.1),
        axis.line = element_line(size = 0.1),
        plot.title = element_text(size = 5)) +
  xlim(c(1.5, 7)) +
  ggtitle('pearson R and P from longest to shortest UTR:
          \n R -0.1506148, P 2.24e-11
          \n R -0.1105765, P 1.25e-06
          \n R -0.1395903, P 1.42e-08
          \n R -0.1232038, P 5.46e-05')

pdf("S6A_are.tpu.hl.pdf", width = 3, height = 3)
are.tpu.hl 
dev.off()

#-------------S6B: ELAVL deseq vs. hl

# Wilcoxon rank-sum test (Mann-Whitney, paired = F)

# Define a vector to hold the p-values
p_values <- c()

# Perform tests and store the p-values
p_values <- c(p_values, (data %>% 
                           filter(percentile_are == '(0,1]') %>%
                           {wilcox.test(hl.GFP ~ elavl, data = ., paired = F)})$p.value) # 25% ARE: p-value = 0.0001891

p_values <- c(p_values, (data %>% 
                           filter(percentile_are == '(1,3]') %>%
                           {wilcox.test(hl.GFP ~ elavl, data = ., paired = F)})$p.value) # 50% ARE: p-value = 0.0008707

p_values <- c(p_values, (data %>% 
                           filter(percentile_are == '(3,8]') %>%
                           {wilcox.test(hl.GFP ~ elavl, data = ., paired = F)})$p.value) # 75% ARE: p-value = 3.624e-06

p_values <- c(p_values, (data %>% 
                           filter(percentile_are == '(8,95]') %>%
                           {wilcox.test(hl.GFP ~ elavl, data = ., paired = F)})$p.value) # 100% ARE: p-value = 0.06219

# Adjust the p-values using the Benjamini-Hochberg procedure
p_values_adjusted <- p.adjust(p_values, method = "BH")
print(p_values_adjusted)

# plot
elavl.are.plot <- 
  data %>% 
  ggplot(aes(x = hl.GFP,
             y = percentile_are,
             fill = elavl)) +
  geom_boxplot(aes(fill = elavl), lwd = 0.15, outlier.shape = NA) +
  scale_fill_manual(values = c("red", 'grey')) +
  theme_bw() +
  theme(axis.text.x = element_text(angle = 0, vjust = 0.5, hjust=1),
        axis.text.y = element_text(size=5),
        axis.ticks = element_line(size = 0.1),
        axis.line = element_line(size = 0.1),
        plot.title = element_text(size = 5)) +
  xlim(c(1, 9)) +
  ggtitle('Wilcoxon runk-sum pvalues Elavl.down vs. none, \n from highest to lowest ARE percentile:
          \n 0.062
          \n 1.45e-05
          \n 1.16e-03 
          \n 3.78e-04')

# Save the plot to a PDF file
pdf("S6B_elavl.are.pdf", width = 3, height = 3)
elavl.are.plot
dev.off()
