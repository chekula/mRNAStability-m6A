library (tidyverse)
library (dplyr)
library (ggplot2)


# load input hl data
hl.GFP <- read_tsv('../data/hl_table_full_niko.tsv') %>% 
  filter(filtered == 'pseudoR2',
         condition == 'GFP') %>% 
  select(gene_id,
         hl_ftc3_tc,
         condition) %>% 
  rename(hl = hl_ftc3_tc) %>% 
  drop_na()

hl.dnCaf1 <- read_tsv('../data/hl_table_full_niko.tsv') %>% 
  filter(filtered == 'pseudoR2',
         condition == 'dnCaf1') %>% 
  select(gene_id,
         hl_ftc3_tc,
         condition) %>% 
  rename(hl = hl_ftc3_tc) %>% 
  drop_na()

# write function to compute RNA probabily to localize based on hl using Fonkeu et al. formula.
# this function take a dataframe with an 'hl' column as input, 
# and it will return a dataframe with an additional 'Rdens_std' column computed for 
# each combination of 'hl' and 'x' values, where x is neurite distance to which mRNA localizes
# D: diffusion coefficient for mRNA
# V: transport velocity of mRNA granule
# K: degradation rate of mRNA
# b: transcription rate of mRNA
# x: distance to which mRNA localizes

compute_rdens <- function(df) {
  # Constants
  D <- 3.0 * 10^(-3)
  V <- 5.8 * 10^(-3)
  b <- 0.0109
  
  # Create a vector 'x' from 1 to 500 with a step of 5
  x <- seq(1, 500, by = 5)
  
  # Generate all combinations of 'hl' and 'x'
  df_expanded <- merge(df, data.frame(x = x), all = TRUE)
  
  # Calculate K and lambda for each hl value
  df_expanded$K <- log(2) / (df_expanded$hl * 3600)
  df_expanded$lambda <- (sqrt(V^2 + 4 * df_expanded$K * D) - V) / (2 * D)
  
  # Calculate Rdens for each x and hl combination
  df_expanded$Rdens <- (b * df_expanded$lambda / df_expanded$K) * exp(-df_expanded$lambda * df_expanded$x)
  
  # Standardize Rdens by dividing each element by the corresponding Rdens value for x = 1
  df_expanded <- df_expanded %>%
    group_by(hl) %>%
    mutate(Rdens_std = Rdens / Rdens[x == 1]) %>%
    ungroup()
  
  return(df_expanded)
}

# compute RNA probabilities to localize base don hl using compute_rdens
RNAdens.GFP <- compute_rdens(hl.GFP)
RNAdens.dnCaf1 <- compute_rdens(hl.dnCaf1)

write_csv(RNAdens.GFP , 'RNAdens.GFP.csv')
write_csv(RNAdens.dnCaf1 , 'RNAdens.dnCaf1.csv')

# -------------
# tile plot 
colnames(RNAdens.GFP)

# ranking for hl distribution
RNAdens.GFP_ranked <- RNAdens.GFP  %>% 
  drop_na() %>% 
  arrange(hl) %>% 
  mutate(rank = dense_rank(hl))

RNAdens.dnCaf1_ranked <- RNAdens.dnCaf1  %>% 
  drop_na() %>% 
  arrange(hl) %>% 
  mutate(rank = dense_rank(hl))

p1 <- 
  RNAdens.GFP_ranked %>% 
  ggplot(aes(x = x, 
             y = rank)) + 
  geom_tile(aes(fill = Rdens_std)) + 
  scico::scale_fill_scico(palette = "lajolla", begin = 0, name = "Probability:") +
  labs(y = "transcripts arranged by hl, 11042", x = "Neurite distance, um") +
  theme(aspect.ratio = 1, 
        panel.grid = element_blank(), 
        panel.border = element_blank(),
        axis.text.y = element_blank(),
        axis.ticks.y = element_blank())+
  ggtitle('GFP')

p2 <- 
  RNAdens.dnCaf1_ranked %>% 
  ggplot(aes(x = x, 
             y = rank)) + 
  geom_tile(aes(fill = Rdens_std)) + 
  scico::scale_fill_scico(palette = "lajolla", begin = 0, name = "Probability:") +
  labs(y = "transcripts arranged by hl, 10368", x = "Neurite distance, um") +
  theme(aspect.ratio = 1, 
        panel.grid = element_blank(), 
        panel.border = element_blank(),
        axis.text.y = element_blank(),
        axis.ticks.y = element_blank())+
  ggtitle('dnCaf1')

pdf('3C_tile.plots.pdf', width = 3, height = 6)
grid.arrange(p1, p2)
dev.off()



