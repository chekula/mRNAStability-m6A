library (tidyverse)
library (dplyr)
library (ggplot2)
library(gridExtra)
library(rstatix)
library('ggforce')


# load the data

loc <- read.csv('../data/multiple_loc_datasets.csv') 

annotation <- loc %>% 
  select(dataset_loc,
         id) %>% 
  aggregate(. ~ dataset_loc, FUN = unique) %>% 
  arrange(id)

write.csv(annotation, '4D_datasets_annotation.csv', row.names = F)


# ++++++++++++++++ codon optimality across datasets

# histo: tAI counts are normally ditributed -> Welch test
loc %>% 
  ggplot(aes(tAI)) +
  geom_histogram(aes(fill = loc_var), lwd = 0.15, outlier.size = 0.1, outlier.color = 'white') +
  facet_grid(vars(id))


p <- loc %>% 
  ggplot(aes(tAI)) +
  geom_boxplot(aes(fill = loc_var), lwd = 0.15, outlier.size = 0.1, outlier.color = 'white') +
  facet_grid(vars(id)) +
  scale_fill_manual(values = c('forestgreen', 'grey', "royalblue")) +
  xlim(c(0.17, .25)) +
  theme_bw() +
  theme(text = element_text(size = 12),
        axis.text = element_text(color = 'black', size = 12),
        axis.text.x = element_text(size = 5, angle = 0),
        axis.text.y = element_text(size = 5),
        axis.ticks = element_line(size = 0.1),
        axis.line = element_line(size = 0.1),
        plot.title = element_text(size = 5),
        panel.grid.major = element_blank(),
        panel.grid.minor = element_blank())

pdf("6B_datasets.tAI.pdf", width = 3, height = 3)
p
dev.off()


# Welch test

loc %>%
  filter(id == 1, loc_var != 'Not significant') %>%
  with(t.test(tAI ~ loc_var, data = ., paired = FALSE, var.equal = FALSE))

loc  %>%
  filter(id == 2, loc_var != 'Not significant') %>%
  with(t.test(tAI ~ loc_var, data = ., paired = FALSE, var.equal = FALSE))

loc  %>%
  filter(id == 3, loc_var != 'Not significant') %>%
  with(t.test(tAI ~ loc_var, data = ., paired = FALSE, var.equal = FALSE))

loc  %>%
  filter(id == 4, loc_var != 'Not significant') %>%
  with(t.test(tAI ~ loc_var, data = ., paired = FALSE, var.equal = FALSE))

loc %>%
  filter(id == 5, loc_var != 'Not significant') %>%
  with(t.test(tAI ~ loc_var, data = ., paired = FALSE, var.equal = FALSE))

loc  %>%
  filter(id == 6, loc_var != 'Not significant') %>%
  with(t.test(tAI ~ loc_var, data = ., paired = FALSE, var.equal = FALSE))

loc  %>%
  filter(id == 7, loc_var != 'Not significant') %>%
  with(t.test(tAI ~ loc_var, data = ., paired = FALSE, var.equal = FALSE))

loc  %>%
  filter(id == 8, loc_var != 'Not significant') %>%
  with(t.test(tAI  ~ loc_var, data = ., paired = FALSE, var.equal = FALSE))

loc  %>%
  filter(id == 9, loc_var != 'Not significant') %>%
  with(t.test(tAI  ~ loc_var, data = ., paired = FALSE, var.equal = FALSE))

loc  %>%
  filter(id == 10, loc_var != 'Not significant') %>%
  with(t.test(tAI ~ loc_var, data = ., paired = FALSE, var.equal = FALSE))

loc  %>%
  filter(id == 11, loc_var != 'Not significant') %>%
  with(t.test(tAI  ~ loc_var, data = ., paired = FALSE, var.equal = FALSE))

loc  %>%
  filter(id == 12, loc_var != 'Not significant') %>%
  with(t.test(tAI ~ loc_var, data = ., paired = FALSE, var.equal = FALSE))

loc  %>%
  filter(id == 13, loc_var != 'Not significant') %>%
  with(t.test(tAI ~ loc_var, data = ., paired = FALSE, var.equal = FALSE))
