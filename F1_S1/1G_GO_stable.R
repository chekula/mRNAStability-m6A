library (tidyverse)
library (dplyr)
library (ggplot2)
library(gridExtra)
library(rstatix)
library(ggpubr)


## load and merge HL and Deseq data

hl <- read_tsv('../data/hl_table_full_niko.tsv') %>% 
  filter(filtered == 'pseudoR2', condition == 'GFP') %>% 
  rename(hl = hl_ftc3_tc) %>% 
  filter(!is.na(hl)) %>% 
  select(gene_id, hl) 

deseq <- read.csv('../data/SupplementaryTable3_DE.csv') %>% 
  filter(experiment == 'PCN.ThreeWay_dnCaf1.CompSep') %>% 
  filter(comparison == 'Pro_GFP.vs.Cyto_GFP',
         DE == 'Projections') %>%   
  select(gene_id, gene_name, log2FoldChange, DE)

data <- merge(deseq, hl, by = 'gene_id') %>% 
  arrange(desc(hl)) 

data.GO.top10 <- data[1:round(nrow(data)*0.1),]

write.csv(data.GO.top10, 'data.GO.top10.csv')

# use 'data.GO.top10.csv' as input into gProfiler to obtain 'top10_gProfiler_mmusculus_13-03-2023_12-48-33__intersections.csv'

# plot GO terms for neurite-localized top 10% stable mRNAs

go <- read.csv('top10_gProfiler_mmusculus_13-03-2023_12-48-33__intersections.csv')

go.plot <- 
  go %>%
  filter(term_name %in% c('structural constituent of ribosome', 
                          'cellular respiration',
                          'translation', 
                          'ribosome', 
                          'mitochondrial protein-containing complex',
                          'neuron projection', 
                          'neuron to neuron synapse')) %>% 
  ggplot(aes(x = intersection_query_percent, y = reorder(term_name, intersection_query_percent))) +
  geom_col(aes(fill = negative_log10_of_adjusted_p_value), size = 0.25, color = "black") +
  scale_fill_gradientn(colors = RColorBrewer::brewer.pal(9, "Reds")) +
  labs(x = "Intersection / query, %", y = 'Term name') +
  theme_bw() +
  theme(axis.text.x = element_text(angle = 0, vjust = 0.5, hjust = 1),
        text = element_text(size = 5),
        axis.text = element_text(color = 'black', size = 5),
        plot.title = element_text(size = 5)) 

pdf("1G_GO_figure.pdf", width = 4, height = 3)
go.plot 
dev.off()

