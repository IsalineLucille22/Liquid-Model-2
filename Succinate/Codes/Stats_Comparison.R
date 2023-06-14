library(ggplot2)
library(ggpubr)
library(readxl)
library(dplyr)
library(grid)
library(Hmisc)

rm(list = ls())

setwd("/Users/isalinelucille-guex/Documents/Liquid-models/Succinate/Codes/Data")

data = read_excel("DensitiesComparisonManuRenyi.xlsx",sheet = 10)

#Show a random sample (always the same with set.seed function)
set.seed(1234)
dplyr::sample_n(data,10)

#Generate frequency tables
#Violin plot
table(data$Ratio, data$Model)

bp = ggplot(data, aes(x = Model, y = ProportionsPVE, color = Model )) + geom_violin() + geom_point(position = position_jitter(seed = 1, width = 0.2)) + theme(legend.position = "none")
stat_compare_means(method = "t.test", aes(group = Model), label.y = 0.9)
bp = bp + stat_summary(fun = mean, geom = "crossbar", shape = 3, size = 0.1, col = "gray")
bp = bp + stat_summary(fun = median, geom = "crossbar", shape = 3, size = 0.1, col = "black")
bp + scale_color_manual(values = c("dodgerblue4", "dodgerblue", "deepskyblue", "firebrick2"))


#Boxplots

data$Ratio = factor(data$Ratio,levels = c("100:1", "10:1", "1:1", "1:10", "1:100"))

bp = ggboxplot(data, x = "Ratio", y = "ProportionsPVE",
               color = "Model", palette = "jco") + 
  stat_compare_means(method = "t.test", aes(group = Model), label.y = 0.9)

bp + scale_color_manual(values = c("dodgerblue4", "dodgerblue", "deepskyblue", "firebrick2"))

