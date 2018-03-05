# Title     : TODO
# Objective : TODO
# Created by: schencro
# Created on: 3/2/18

library(ggplot2)

setwd("/Users/schencro/Desktop/Oxford/Rotation_1/CNN/Model/ModelComparisons/")

df <- read.csv("CellAucComparisons.csv", sep=",", stringsAsFactors = F, header = T)
df$Dataset <- as.factor(df$Dataset)
fit <- aov(AUC ~ Dataset, df)
summary(fit)
TukeyHSD(fit)

p <- ggplot() + geom_boxplot(data=df, aes(x=Dataset, y=AUC, color=Dataset), notch=TRUE) + 
  geom_point(data=df, aes(x=Dataset, y=AUC),color="black", alpha=0.2, inherit.aes = F) +
  geom_line(data=df, aes(x=Dataset, y=AUC, group=NormCell), alpha=0.1) + 
  theme_minimal() + guides(fill=FALSE, color=FALSE) + xlab("") + scale_color_brewer(palette="Set2") +
  ggtitle("Individual AUC Differences")
p
