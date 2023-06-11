library (tidyverse)
library (dplyr)
library (ggplot2)
library(gridExtra)
library(rstatix)

## load and merge HL and Deseq data

hl <- read_tsv('../data/hl_table_full_niko.tsv') %>% 
  filter(filtered == 'pseudoR2') %>% 
  select(gene_id,
         hl_ftc3_tc,
         condition) %>% 
  pivot_wider(values_from = hl_ftc3_tc,
              names_from = condition, names_prefix = "hl.") %>% 
  mutate(delta.hl = hl.dnCaf1 - hl.GFP)

deseq <- read.csv('../data/SupplementaryTable3_DE.csv') %>% 
  filter(comparison == "Pro_GFP.vs.Cyto_GFP" |
           comparison == "Pro_dnCaf1.vs.Cyto_dnCaf1"  )

merge <- merge(deseq, hl, by = 'gene_id') %>% 
  select(gene_id,
         log2FoldChange,
         comparison,
         hl.dnCaf1,
         hl.GFP,
         delta.hl) %>% 
  pivot_wider(names_from = comparison,
              values_from = log2FoldChange) %>% 
  drop_na() %>% 
  mutate(delta.loc = Pro_dnCaf1.vs.Cyto_dnCaf1 - Pro_GFP.vs.Cyto_GFP)
  
  # plotting percentiles of delta HL/hl.gfp  against delta localization for hl.dnCaf1 > 5.6 hr
  
  # Pearson cor: R 0.1335041, p-value = 3.764e-14
  
  merge.more5.6 <- merge %>% 
    filter(hl.dnCaf1 > 5.6) %>% 
    mutate(delta.hl.perc = delta.hl/hl.GFP*100) %>% 
    mutate(percentile.delta.hl.perc = cut(delta.hl.perc, 
                                     quantile(delta.hl.perc, seq(0,1, .2), 
                                              na.rm = T), 
                                     include_lowest = TRUE))
  
  delta.hl <- merge.more5.6 %>% 
    select(delta.hl) %>% 
    unlist()
  
  delta.loc <- merge.more5.6 %>% 
    select(delta.loc) %>% 
    unlist()
  
  Pearson <- cor.test(delta.hl, delta.loc, method = 'pearson')
  
  # plotting
  caf.hl.perc.more5.6.plot <- 
    merge.more5.6 %>% 
    drop_na() %>% 
    ggplot(aes(percentile.delta.hl.perc,
               delta.loc)) +
    geom_boxplot(aes(fill = percentile.delta.hl.perc), lwd=0.25, outlier.size = 0.1, outlier.shape =  NA) +
    scale_fill_manual(values = alpha(rep("red", 5), seq(0.2, 0.9, length.out = 5))) +
    theme_bw() +
    theme(axis.text.x = element_text(angle = 90, vjust = 0.5, hjust=1),
          axis.text.y = element_text(size=5),
          axis.ticks = element_line(size = 0.1),
          axis.line = element_line(size = 0.1),
          plot.title = element_text(size = 5)) +
    ylim(c(-2, 2)) +
    ggtitle('hl.dncaf1 > 5.6 hr,
    \n Pearson R 0.1335041 
          \n p = 3.764e-14')+
    labs(x = 'deltaHL, %')
  
  pdf("3E_d.loc.vs.d.hl.pdf", width = 3, height = 3)
  caf.hl.perc.more5.6.plot
  dev.off()
  
  # plotting percentiles of delta HL against delta localization for hl.dnCaf1 < 5.6 hr
  
  # ---------- hl.dnCaf1 < 5.6 hr
  # Pearson cor: 0.1061226, p-value = 6.404e-11
  merge.less5.6 <- merge %>% 
    filter(hl.dnCaf1 < 5.6) %>% 
    mutate(delta.hl.perc = delta.hl/hl.GFP*100) %>% 
    mutate(percentile.delta.hl.perc = cut(delta.hl.perc, 
                                     quantile(delta.hl.perc, seq(0,1, .2), 
                                              na.rm = T), 
                                     include_lowest = TRUE))
  
  delta.hl <- merge.less5.6  %>% 
    select(delta.hl) %>% 
    unlist()
  
  delta.loc <- merge.less5.6 %>% 
    select(delta.loc) %>% 
    unlist()
  
  Pearson <- cor.test(delta.hl, delta.loc, method = 'pearson')
  
  
  caf.hl.perc.less5.6.plot <- 
    merge.less5.6 %>% 
    filter(hl.dnCaf1 < 5.6) %>% 
    drop_na() %>% 
    ggplot(aes(percentile.delta.hl.perc,
               delta.loc)) +
    geom_boxplot(aes(fill = percentile.delta.hl.perc), lwd=0.25, outlier.size = 0.1, outlier.shape =  NA) +
    scale_fill_manual(values = alpha(rep("red", 5), seq(0.2, 0.9, length.out = 5))) +
    theme_bw() +
    theme(axis.text.x = element_text(angle = 90, vjust = 0.5, hjust=1),
          axis.text.y = element_text(size=5),
          axis.ticks = element_line(size = 0.1),
          axis.line = element_line(size = 0.1),
          plot.title = element_text(size = 5)) +
    ylim(c(-2, 2)) +
    ggtitle('hl.dnCaf1 < 5.6 hr,
            \n Pearson R 0.1061226  
            \n p = 6.404e-11')+
    labs(x = 'deltaHL, %')
  
  pdf("S2B_d.loc.vs.d.hl.pdf", width = 3, height = 3)
  caf.hl.perc.less5.6.plot
  dev.off()
  
  