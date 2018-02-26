# Title     : ROC Curves
# Objective : Plot ROC and AUCs for model testing data.
# Created by: schencro
# Created on: 2/25/18

# required
require(ggplot2, quietly = TRUE)
require(RColorBrewer, quietly = TRUE)

args <- commandArgs(trailingOnly = TRUE)

# Parse arguments
input.dir <- args[1]  # Input Directory
out.dir <- args[2]  # Output Directory
plot.prefix <- args[3]  # prefix

print(input.dir)
print(out.dir)
print(plot.prefix)

inputFile <- paste(input.dir, "/roc_curve_data.csv", sep="")
df <- read.csv(inputFile, header=T, sep=',')
inputFile <- paste(input.dir, "/", plot.prefix, '.roc_curve_data.micromacro.csv', sep="")
dfAves <- read.csv(inputFile, header=T, sep=",")

min <- paste("Min. AUC in Cell ", unique(subset(df,df$AUC==summary(df$AUC)[1])$Cell),": ", round(summary(df$AUC)[1], 4), sep="")
max <- paste("Max. AUC in Cell ", unique(subset(df,df$AUC==summary(df$AUC)[6])$Cell),": ", round(summary(df$AUC)[6], 4), sep="")
macro <- paste("Macro-Average (AUC: ",round(unique(subset(dfAves, dfAves$SumType=="macroall")$AUC),4),")", sep ="")
micro <- paste("Micro-Average (AUC: ",round(unique(subset(dfAves, dfAves$SumType=="microall")$AUC),4),")", sep ="")

# Overall Macro Micro
P <- ggplot() + geom_line(data=df, aes(x=fpr, y=tpr, group=Cell), alpha=0.2, size=0.2, color="grey", inherit.aes = F) +
  geom_line(data=dfAves, aes(x=fpr, y=tpr, colour=SumType), size=0.5, inherit.aes = F) +
  theme_minimal() +
  scale_color_manual(labels = c(macro, micro), values = c("darkorange2", "chartreuse3")) +
  xlab("False Positive Rate (1-Specificity)") + ylab("True Positive Rate (Sensitivity)") +
  ggtitle(paste(min, "  ", max, sep="")) + geom_abline(slope=1, intercept=0, size=0.2, linetype="dashed") +
  theme(legend.position = c(0.8, 0.25)) + guides(color=guide_legend(title=""))
ggsave(P, filename = paste(out.dir, plot.prefix, ".ROCCurve.png", sep=""), width = 10, height = 10)