library(ggplot2)
library(ggpubr)
library(readxl)
library(dplyr)
library(grid)
library(Hmisc)

rm(list = ls())

# setwd("/Users/isalinelucille-guex/switchdrive/shared UNIL/Ppu_Pve_growth_data/version 2/figures/Overleaf/Fig6_BoxPlotComp")
setwd("/Users/isalinelucille-guex/switchdrive/shared UNIL/Ppu_Pve_growth_data/version 2/figures/Overleaf/Fig6_EffectActivationFunction/BoxPlotComp")

data = read_excel("DensitiesComparisonManuRenyi.xlsx",sheet = 10)

data = data[which(data$Model != "M3"), ]
data = data[which(data$Model != "M1"), ]

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

# bp = ggboxplot(data_A, x = "Ratio", y = "ProportionsPVE",
#                color = "Model", palette = "jco") +
#   stat_compare_means(method = "t.test", aes(group = Model), label.y = 0.9) + geom_violin(trim=FALSE)

#Stationary Pve proportions
data$Ratio = factor(data$Ratio,levels = c("100:1", "10:1", "1:1", "1:10", "1:100"))

bp = ggboxplot(data, x = "Ratio", y = "ProportionsPVE",
               color = "Model", palette = "jco") + 
  stat_compare_means(method = "t.test", aes(group = Model), label.y = 0.9)

bp + scale_color_manual(values = c("dodgerblue4", "dodgerblue", "deepskyblue", "firebrick2"))

#Sum stationary absolute abundances
bp = ggboxplot(data, x = "Ratio", y = "SumAbundances",
               color = "Model", palette = "jco") + 
  stat_compare_means(method = "t.test", aes(group = Model), label.y = 1e-04)

bp + scale_color_manual(values = c("dodgerblue4", "dodgerblue", "deepskyblue", "firebrick2"))


############Boxplots for supplementary information

rm(list = ls())

setwd("/Users/isalinelucille-guex/switchdrive/shared UNIL/Ppu_Pve_growth_data/version 2/figures/Overleaf/Fig6_EffectActivationFunction/BoxPlotComp")

data = read_excel("DensitiesComparisonManuRenyi.xlsx",sheet = 11)

data = data[which(data$Model != "S1"), ]

#Show a random sample (always the same with set.seed function)
set.seed(1234)
dplyr::sample_n(data,10)

#Generate frequency tables
#Violin plot
table(data$Ratio, data$Model)


data$Ratio = factor(data$Ratio,levels = c("100:1", "10:1", "1:1", "1:10", "1:100"))

# #"ProportionsPVE"
# bp = ggboxplot(data, x = "Ratio", y = "ProportionsPVE",
#                color = "Model", palette = "jco") + 
#   stat_compare_means(method = "t.test", aes(group = Model), label.y = 0.9)

bp = ggboxplot(data, x = "Ratio", y = "SumAbundances",
               color = "Model", palette = "jco") + 
  stat_compare_means(method = "t.test", aes(group = Model), label.y = 8e-05)

bp + scale_color_manual(values = c("firebrick2", "dodgerblue4", "dodgerblue"))



#Two-way ANOVA

res.anova = aov(Proportions ~ Ratio + Model, data = data_A)
summary(res.anova)


#Comparison between T-test with two level (Obs or M) and One-way Anova. 

data_Tt = data_A[which(data_A$Ratio == "1:100"), ]
x = as.numeric(unlist(data_Tt[which(data_Tt$Model == "M1"), 2]))
y = as.numeric(unlist(data_Tt[which(data_Tt$Model == "Obs"), 2]))
ks.test(x, y)

#Shapiro-Wilk test to check the normality
#The Shapiro–Wilk test tests the null hypothesis that a sample x1, ..., xn came from a normally distributed population. The test statistic is 

with(data_Tt, shapiro.test(Proportions[Model == "M1"]))
with(data_Tt, shapiro.test(Proportions[Model == "Obs"]))

#Bartlett's test in order to test the homogeneity of variances. Assumption of normality must be satisfied.

res.Bart = bartlett.test(Proportions ~ Model, data=data_Tt)
res.Bart

res.Var = var.test(Proportions ~ Model, data=data_Tt)
res.Var

#Perform t-test
res.Ttest = t.test(Proportions ~ Model, data = data_Tt)
res.Ttest

Obs_1 = as.numeric(unlist(data_Tt[which(data_Tt$Model == "Obs"), 2]))
Obs_2 = as.numeric(unlist(data_Tt[which(data_Tt$Model == "M2"), 2]))
X_1 = mean(Obs_1)
X_2 = mean(Obs_2)
Var_1 = var(Obs_1)
Var_2 = var(Obs_2)

t = (X_1-X_2)/sqrt((Var_1+Var_2)/4)
t.test(Obs_1,Obs_2)

#Peform one-way anova with 3 levels of "Model" ("Obs", "M1", "M2")

data_aov_3 = data[which(data$Ratio == "100:1"), ]

#Shapiro-Wilk test to check the normality

with(data_aov_3, shapiro.test(Proportions[Model == "Obs"]))
with(data_aov_3, shapiro.test(Proportions[Model == "M1"]))
with(data_aov_3, shapiro.test(Proportions[Model == "M2"]))

#Bartlett's test in order to test the homogeneity of variances. Assumption of normality must be satisfied.

res.Bart = bartlett.test(Proportions ~ Model, data=data_aov_3)
res.Bart

#Perform Anova 

res.aov_3= aov(Proportions ~ Model, data = data_aov_3)
summary(res.aov_3)

#Perform LSD test

M_Obs = mean(unlist(data_aov_3[which(data_aov_3$Model == "Obs"), 2]))
M_M1 = mean(unlist(data_aov_3[which(data_aov_3$Model == "M1"), 2]))
M_M2 = mean(unlist(data_aov_3[which(data_aov_3$Model == "M2"), 2]))

LSD = qt(1-0.05/2, 9)*sqrt(0.00042*(1/4 + 1/4))

diff_ObsM1 = abs(M_Obs - M_M1) - LSD
diff_ObsM2 = abs(M_Obs - M_M2) - LSD

#Example found in R-manual

x <- rnorm(50)
y <- runif(30)
# Do x and y come from the same distribution?
ks.test(x, y)
# Does x come from a shifted gamma distribution with shape 3 and rate 2?
ks.test(x+2, "pgamma", 3, 2) # two-sided, exact
ks.test(x+2, "pgamma", 3, 2, exact = FALSE)
ks.test(x+2, "pgamma", 3, 2, alternative = "gr")

# test if x is stochastically larger than x2
x2 <- rnorm(50)
plot(ecdf(x), xlim = range(c(x, x2)))
plot(ecdf(x2), add = TRUE, lty = "dashed")
t.test(x, x2, alternative = "g")
wilcox.test(x, x2, alternative = "g")
ks.test(x, x2, alternative = "l")

