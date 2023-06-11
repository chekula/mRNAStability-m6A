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


# ++++++++++++++++ m6A across datasets

# calculate % of m6_enriched in N-, S- and not significant transcripts
data.m6a <- loc %>% 
  group_by(id, loc_var) %>% 
  summarize(m6a.enrich.percentage = sum(m6a.cons.enrich == 'm6A enriched') / n()*100)

p <- 
  data.m6a %>%
  ggplot(aes(x = id, y = m6a.enrich.percentage, fill = loc_var)) +
  geom_col(position = 'dodge') +
  scale_fill_manual(values = c('forestgreen', 'grey', "royalblue")) +
  theme_bw() +
  theme(text = element_text(size = 4),
    axis.text = element_text(color = 'black', size = 4),
    axis.text.x = element_text(size = 4),
    axis.text.y = element_text(size = 4),
    axis.line = element_line(size = 0.1),
    axis.ticks.x = element_line(size = 0.1),
    axis.ticks.y = element_line(size = 0.1),
    panel.grid.major = element_blank(),
    panel.grid.minor = element_blank(),
    plot.title = element_text(size = 4)) +
  scale_x_discrete(limits = unique(data.m6a$id))

  
pdf("4D_datasets.m6a.pdf", width = 3, height = 3)
p
dev.off()


# Chi-squared test to assess frequencies 

loc %>%
  filter(id == 1, loc_var != 'Not significant') %>%
  with(chisq.test(m6a.cons.enrich, loc_var))

loc %>%
  filter(id == 2, loc_var != 'Not significant') %>%
  with(chisq.test(m6a.cons.enrich, loc_var))

loc %>%
  filter(id == 3, loc_var != 'Not significant') %>%
  with(chisq.test(m6a.cons.enrich, loc_var))

loc %>%
  filter(id == 4, loc_var != 'Not significant') %>%
  with(chisq.test(m6a.cons.enrich, loc_var))

loc %>%
  filter(id == 5, loc_var != 'Not significant') %>%
  with(chisq.test(m6a.cons.enrich, loc_var))

loc %>%
  filter(id == 6, loc_var != 'Not significant') %>%
  with(chisq.test(m6a.cons.enrich, loc_var))

loc %>%
  filter(id == 7, loc_var != 'Not significant') %>%
  with(chisq.test(m6a.cons.enrich, loc_var))

loc %>%
  filter(id == 8, loc_var != 'Not significant') %>%
  with(chisq.test(m6a.cons.enrich, loc_var))

loc %>%
  filter(id == 9, loc_var != 'Not significant') %>%
  with(chisq.test(m6a.cons.enrich, loc_var))

loc %>%
  filter(id == 10, loc_var != 'Not significant') %>%
  with(chisq.test(m6a.cons.enrich, loc_var))

loc %>%
  filter(id == 11, loc_var != 'Not significant') %>%
  with(chisq.test(m6a.cons.enrich, loc_var))

loc %>%
  filter(id == 12, loc_var != 'Not significant') %>%
  with(chisq.test(m6a.cons.enrich, loc_var))

loc %>%
  filter(id == 13, loc_var != 'Not significant') %>%
  with(chisq.test(m6a.cons.enrich, loc_var))
