library (tidyverse)
library (dplyr)
library (ggplot2)
library(gridExtra)
library(rstatix)
library('ggforce')


# load predictions
predictions_output <- read.csv('7C_loc.ranger_predictions.csv')

# vizualize distribution of actual and predicted localisation (log2FC)

predictions_output %>% 
  ggplot() +
  geom_histogram(aes(ranger_p.prediction, fill = 'red'), alpha = 0.5) +
  geom_histogram(aes(log2FoldChange, fill = log2FoldChange), alpha = 0.5) 

# select transcripts predicted to be localizaed (log2FC > 0.58) and save for GOterm analysis with Gprofiler

predicted.neu.loc <- 
  predictions_output %>% 
  filter(ranger_p.prediction > 0.58)

write.csv(predicted.neu.loc, "all.predicted.neu.loc.csv")

# plotting selected GO terms for transcripts predicted to be neurite-localized by the ranger model.

go <- read.csv('../data/Predicted.neu_gProfiler_mmusculus_16-03-2023_10-53-25__intersections.csv')

go.plot <- 
  go %>%
  filter(term_name %in% c('structural constituent of ribosome', 
                          'large ribosomal subunit',
                          'cytoplasmic translation', 
                          'cytosolic ribosome', 
                          'inner mitochondrial membrane protein complex',
                          'synapse')) %>% 
  ggplot(aes(x = intersection_query_percent, y = reorder(term_name, intersection_query_percent))) +
  geom_col(aes(fill = negative_log10_of_adjusted_p_value), size = 0.25, color = "black") +
  scale_fill_gradientn(colors = RColorBrewer::brewer.pal(9, "Reds")) +
  labs(x = "Intersection / query, %", y = 'Term name') +
  theme_bw() +
  theme(axis.text.x = element_text(angle = 0, vjust = 0.5, hjust = 1),
        text = element_text(size = 5),
        axis.text = element_text(color = 'black', size = 5),
        plot.title = element_text(size = 5)) +
  ggtitle('top6 GOterms by p-value')

pdf("7D_GOterms.pdf", width = 4, height = 2)
go.plot 
dev.off()

# percentgae of correctly predicted N-localized transcripts

subset <- subset(predictions_output, predictions_output$log2FoldChange > 0.58)
nrow(subset)

subset2 <- subset(predictions_output, predictions_output$log2FoldChange > 0.58 &
                    predictions_output$ranger_p.prediction > 0.58)
nrow(subset2)

correctly.pred = 142/533*100
