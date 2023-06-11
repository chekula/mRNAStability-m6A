library (tidyverse)
library (dplyr)
library (ggplot2)
library(gridExtra)


## load HL 

hl <- read.csv('../data/HL.csv') %>% 
  filter(group %in% c("Pro","Cyto", 'Neurite', "Soma")) %>% 
  select(gene_id,
        Half_life,
        group) %>% 
  pivot_wider(names_from = group,
              values_from = Half_life,
              names_prefix = 'hl.') %>% 
  mutate(delta.hl = hl.Pro - hl.Cyto,
         delta.hl.mESC = hl.Neurite - hl.Soma) 

hl %>% 
  filter(!is.na(delta.hl)) %>% 
  summarize(median.delta.hl = median(delta.hl),
            mean.delta.hl = mean(delta.hl))

hl %>%
  filter(!is.na(delta.hl)) %>%
  summarize(median = median(delta.hl),
            mean = mean(delta.hl))

#  histo 
histo <- hl %>% 
  ggplot(aes(x = delta.hl)) +
  geom_histogram(bins = 50,
                 color = 'black',
                 lwd = .2,
                 alpha = 0.6,
                 position="identity") +
  scale_fill_manual(values = c('red', 'darkgrey')) +
  theme(text = element_text(size=6), 
        axis.text	=element_text(size=6), 
        plot.title = element_text(size=6),
        axis.ticks = element_line(linewidth = 0.15), 
        axis.line = element_line(linewidth = 0.15)) +
  theme_bw() +
  geom_vline(xintercept = 1.06, color = 'red', linetype = 'dashed') +
  xlim(c(-6, 10))

# save as pdf
pdf("1C_delta.hl_histo.pdf", width = 3, height = 3)
histo
dev.off()

