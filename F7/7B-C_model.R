library (tidyverse)
library (dplyr)
library (ggplot2)
library(gridExtra)
library(rstatix)
library('ggforce')
#install.packages("caret", dependencies = c("Depends", "Suggests"))
library(caret)
library(ranger)

# load data (3C PCN, dataset 3)
short.data <- readRDS("../data/short.data-for.modeling.rds")

set.seed(42)

# Shuffle row indices: rows
rows <- sample(nrow(short.data))

# Randomly order data
short.data <- short.data[rows,]

# Determine row to split on: split
split <- round(nrow(short.data)*0.7)

# Create train
train_data <- short.data[1:split,]
saveRDS(train_data, "train.rds")

# Create test
test_data <- short.data[(split+1):nrow(short.data),]
saveRDS(test_data, "test.rds")

train <- train_data %>% 
  select(hl, m6a.enrich, TPU.length, ARE.cts, tAI, log2FoldChange)

test <- test_data %>% 
  select(hl, m6a.enrich, TPU.length, ARE.cts, tAI, log2FoldChange)


# model to predict localization

set.seed(42)

ranger_model <- ranger(log2FoldChange ~ ., data = train, importance = "impurity")

# extract feature importance
importance(ranger_model)

saveRDS(ranger_model, "7C_loc.ranger_model.rds")

ranger_p <- predict(ranger_model, test)

predictions_output <- tibble(ranger_p$prediction, test_data) 

write.csv(predictions_output, 
          "7C_loc.ranger_predictions.csv",
          row.names = F)

predictions_output <- read.csv("7C_loc.ranger_predictions.csv")


rsquared <- R2(predictions_output$ranger_p.prediction, 
                predictions_output$log2FoldChange)

cor <- cor.test(predictions_output$ranger_p.prediction, 
         predictions_output$log2FoldChange, 
         method = 'pearson')

# scatterplot predicted vs. actual localization
p1 <- 
  predictions_output %>% 
  ggplot(aes(x = ranger_p.prediction,
             y = log2FoldChange))+
  geom_hline(yintercept = 0) +
  geom_vline(xintercept = 0) +
  geom_point(alpha = 0.5) +
  geom_smooth()+
  ggtitle('Pearson 0.42, Rsquared 0.18')+
  theme_bw()+
  theme(axis.text.x = element_text(angle = 0, vjust = 0.5, hjust=1))+
  theme(text = element_text(size=8), 
        axis.text	=element_text(color='black', size=8))+
ylim(c(-3,4)) +
  xlim(c(-1,3)) 

pdf("7C_loc.predictions.scatterplot.pdf", width = 3, height = 3)
p1
dev.off()


#   model to predict HL

set.seed(42)

ranger_model <- ranger(hl ~ m6a.enrich+TPU.length+ARE.cts+tAI, data = train, importance = "impurity")

# extract feature importance
importance(ranger_model)

saveRDS(ranger_model, "7B_hl.ranger_model.rds")

ranger_p <- predict(ranger_model, test)

predictions_output <- tibble(ranger_p$prediction, test_data) 

write.csv(predictions_output, 
          "7B_hl.ranger_predictions.csv",
          row.names = F)

predictions_output <- read.csv("7B_hl.ranger_predictions.csv")


rsquared <- R2(predictions_output$ranger_p.prediction, 
              predictions_output$hl)

cor <- cor.test(predictions_output$ranger_p.prediction, 
                predictions_output$hl, 
                method = 'pearson')

# scatter plot predicted vs. actual hl
p2 <- 
  predictions_output %>% 
  ggplot(aes(x = ranger_p.prediction,
             y = hl))+
  geom_hline(yintercept = 0) +
  geom_vline(xintercept = 0) +
  geom_point(alpha = 0.5) +
  geom_smooth()+
  ggtitle('Pearson 0.48, Rsquared 0.23')+
  theme_bw()+
  theme(axis.text.x = element_text(angle = 0, vjust = 0.5, hjust=1))+
  theme(text = element_text(size=8), 
        axis.text	=element_text(color='black', size=8))+
  ylim(c(2,11)) +
  xlim(c(2,11)) 

pdf("7B_hl.predictions.scatterplot.pdf", width = 3, height = 3)
p2
dev.off()


# ---- plotting Gini feature importance

# load the data
importance <- read.csv('../data/feature_importance_loc.csv') 

p3 <- 
  importance  %>% 
  ggplot(aes(x = importance,
             y = loc.predictor))+
  geom_bar(stat = "identity", width = 0.5, color = 'grey') +
  geom_hline(yintercept = 0) +
  geom_vline(xintercept = 0) +
  theme_bw()+
  theme(axis.text.x = element_text(angle = 0, vjust = 0.5, hjust=1))+
  theme(text = element_text(size=8), 
        axis.text	=element_text(color='black', size=8))


importance.hl <- read.csv('../data/feature_importance_hl.csv') 

p4 <- 
  importance.hl  %>% 
  ggplot(aes(x = importance,
             y = hl.predictor))+
  geom_bar(stat = "identity", width = 0.5, color = 'grey') +
  geom_hline(yintercept = 0) +
  geom_vline(xintercept = 0) +
  theme_bw()+
  theme(axis.text.x = element_text(angle = 0, vjust = 0.5, hjust=1))+
  theme(text = element_text(size=8), 
        axis.text	=element_text(color='black', size=8))


pdf("7B-C_feature_importance.pdf", width = 4, height = 4)
grid.arrange(p3, p4, ncol = 1)
dev.off()

