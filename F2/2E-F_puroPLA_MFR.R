library (tidyverse)
library (dplyr)
library (ggplot2)
library(gridExtra)
library(rstatix)
library('ggforce')
library(tidyr)

# 2E: puro-PLA
## load data

puropla <- read.csv('../data/puroPLA.csv')

count_shLarp <- sum(puropla$shRNA == 'shLarp')
count_ctrl <- sum(puropla$shRNA == 'ctrl')


## t-test

statistic <- 
  t.test(puropla$smFISH_per_area ~ puropla$shRNA, paired=F, var.equal = TRUE) 

## violine plots

p1 <- 
  puropla %>% 
  ggplot(aes(x= shRNA,
             y = smFISH_per_area))+
  geom_boxplot(aes(fill = shRNA), alpha = 0.1, color = 'grey')+
  geom_violin(aes(fill = shRNA), alpha = 0.1, color = 'grey')+
  geom_sina(aes(colour = shRNA), size = 2)+
  theme_bw()+
  xlab("shRNA")+
  ylab('puro_PLA per area')+
  scale_color_manual(values = c('darkgrey', 'red'))+
  scale_fill_manual(values = c('darkgrey', 'red'))+
  ggtitle(paste('p-value', statistic[3]))+
  theme(plot.title = element_text(size = 5, face = "bold"))

## saving as pdf
pdf("2E_puropla.pdf", width = 3, height = 3)
p1
dev.off()

# --------------  2F: Weighted mean firing rate (MFR)

#  Load readxl library and WMFR data
library(readxl)

mfr <- read_excel('../data/20230117_Scr_Larp_Elavls_24well_all.xlsx', 
                  sheet = "summary") 

pivoted_mfr <- mfr %>%
  filter(Well != "Treatment") %>%
  pivot_longer(cols = -DIV, 
               names_to = "Well", 
               values_to = "MFR") 

shRNA <- mfr %>% 
  select(-DIV) %>% 
  filter(Well == "Treatment") %>%
  pivot_longer(cols = everything(), 
               names_to = "Well", 
               values_to = "shRNA")
  
combined_mfr <- merge(pivoted_mfr, shRNA, by = "Well") %>% 
  filter(Well != "Well") %>% 
  mutate(MFR = as.numeric(MFR))

write.csv(combined_mfr, '../data/WMFR.csv')

# Permutation test

data <- 
  combined_mfr %>% 
  filter (shRNA != 'Elavls')

# Permutaiton test
# To analyze data over time, a permutation test was performed in R. To generate each permutation dataset, the
# labels of each well (control vs experimental) were randomly shuffled 10,000 times, but data from each well
# were unchanged so as to preserve the correlations between the time points within wells while breaking
# relationships between the group (control vs experimental) and subsequent outcome. A MWU test was
# performed for each of the 10,000 data sets comparing distributions in outcome between groups and a p-value
# was recorded. A permutation p-value for the MWU test was computed as the proportion of permuted data
# MWU p-values that were less than or equal to the MWU p-value from the original un-permuted dataset.

# Function to perform a single permutation
perform_permutation <- function(data) {
  # Shuffle labels
  shuffled_labels <- sample(data$shRNA)
  
  # Assign shuffled labels to each well
  df_permuted <- data %>%
    mutate(shRNA = shuffled_labels)
  
  # Perform MWU test on permuted data
  mwu_test <- wilcox.test(MFR ~ shRNA, data = df_permuted)
  
  return(mwu_test$p.value)
}

# Perform MWU test on original data
mwu_original <- wilcox.test(MFR ~ shRNA, data = data)

# Perform 5,000 permutations
n_permutations <- 5000
permuted_p_values <- replicate(n_permutations, perform_permutation(data))

# Compute permutation p-value
# permuted_p_values <= mwu_original$p.value creates a logical vector with the same length as permuted_p_values, where each element is TRUE if the corresponding permuted p-value is less than or equal to the p-value from the original data, and FALSE otherwise.
# mean(permuted_p_values <= mwu_original$p.value) calculates the proportion of TRUE values in the logical vector, which represents the proportion of permuted datasets with p-values less than or equal to the p-value from the original dataset. This proportion serves as the permutation p-value, which is an estimate of the probability of observing a p-value as extreme or more extreme than the one obtained from the original data, under the null hypothesis (i.e., assuming there is no significant difference between the groups).

permutation_p_value <- mean(permuted_p_values <= mwu_original$p.value)

# Print results
cat("MWU p-value (original data):", mwu_original$p.value, "\n")
cat("Permutation p-value:", permutation_p_value, "\n")


# plotting weighted mean firing rate for Larp KD and Scrambled
p2 <- 
  combined_mfr %>% 
  filter (shRNA != 'Elavls') %>% 
  ggplot(aes(x = interaction(shRNA, DIV),
             y = MFR,
             fill = shRNA)) + 
  geom_boxplot(aes(fill = shRNA), lwd = 0.15, alpha = 0.6)  +
  scale_fill_manual(values = c("red", "grey")) +
  theme_bw()+
  theme(axis.text.x = element_text(angle = 90, vjust = 0.5, hjust=1),
        axis.text.y = element_text(size=5),
        axis.ticks = element_line(size = 0.1),
        axis.line = element_line(size = 0.1),
        plot.title = element_text(size = 5)) +
  ggtitle('Permutation P = 0.016') 

# to show individual datapoints
combined_mfr %>% 
  filter (shRNA != 'Elavls') %>% 
  ggplot(aes(x = interaction(shRNA, DIV),
             y = MFR,
             fill = shRNA)) + 
  geom_boxplot(aes(fill = shRNA), lwd = 0.15, alpha = 0.6, , outlier.shape = NA)  +
  geom_point(aes(color = shRNA), position = position_jitterdodge(), size = 1) +
  scale_fill_manual(values = c("red", "grey")) +
  scale_color_manual(values = c("red", "grey")) +
  theme_bw()+
  theme(axis.text.x = element_text(angle = 90, vjust = 0.5, hjust=1),
        axis.text.y = element_text(size=5),
        axis.ticks = element_line(size = 0.1),
        axis.line = element_line(size = 0.1),
        plot.title = element_text(size = 5)) +
  ggtitle('Permutation P = 0.016') 

pdf("2F_mfr.larp.pdf", width = 3, height = 3)
p2
dev.off()

