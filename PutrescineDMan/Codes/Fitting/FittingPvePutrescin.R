library(ggplot2)
library(ggpubr)
library(readxl)
library(dplyr)
library(grid)
library(Hmisc)
library(FME)
library(deSolve)

## =======================================================================
## Data extraction
## =======================================================================
rm(list = ls())
#dev.off()

setwd("/Users/isalinelucille-guex/Documents/Liquid-models/PutrescineDMan/Codes/Data")

Data_1 = read_excel("DataMaximev4.xlsx", sheet = 7)

TT <- Data_1$...1
Times <- rep(TT, 5)
Data <- data.frame(
  time = Times,
  x = c(Data_1$B6, Data_1$C6, Data_1$D6, Data_1$E6, Data_1$G6)
)

R_0 = 2*2.4*10^(-4)
state = c(x = Data[1,2], y = 0, c = 0, f = 2*2.4*10^(-4), z = Data[1,2]) #Initial conditions used for fitting. These values are drawn from experimental data.

plot(Times, Data[,2])
## =======================================================================
## Model definition
## =======================================================================

#Function used to make the fitting
#The "Monod" model
Monod_2 <- function(t, state, parms){
  with(as.list(c(state, parms)), {
    dx = -kappa_1*x*f + (2*kappa_2 + kappa_3)*y
    dy = kappa_1*x*f - (kappa_2 + kappa_3)*y
    dc = kappa_3*y
    df = -kappa_1*x*f
    list(c(dx, dy, dc, df))
  }) 
}

Monod <- function(t, state, parms){
  with(as.list(c(state, parms)), {
    dx = -kappa_1*x*f + (2*kappa_2 + kappa_3)*y
    dy = kappa_1*x*f - (kappa_2 + kappa_3)*y
    dc = kappa_3*y
    df = -kappa_1*x*f
    dz = dx + dy
    list(c(dx, dy, dc, df, dz))
  }) 
}


##===================================
## Fitted with the "Monod" model #
##===================================
## numeric solution 
## ODEs system
param_init = c(0.5, 0.5, 2*10^5) #Initial guess for the fitting
kappa_2 = param_init[1]
kappa_3 = kappa_2/param_init[2] - param_init[1]
kappa_1 = (kappa_2 + kappa_3)/(0.1*2.4*10^(-4))
parameters_init <- c(kappa_1 = as.numeric(kappa_1), kappa_2 = as.numeric(kappa_2), kappa_3 = as.numeric(kappa_3))

## model cost,
ModelCost <- function(P) {
  out <- ode(y = state, func = Monod, parms = P, times = TT, atol = 1e-10, rtol = 1e-9)
  model = out
  model[,2] = model[,6]
  model = model[,1:2]
  return(modCost(out, Data)) # object of class modCost
}

Fit <- modFit(f = ModelCost, p = parameters_init, lower = rep(0, 3),
              upper = c(5*10^(5), 10, 10)) 

out <- ode(y = state, func = Monod, parms = Fit$par,
           times = TT, atol = 1e-10, rtol = 1e-9)

lines(TT, out[,2], col = "blue", lty = 2)
# legend("right", c("data", "Monod"),
#        lty = c(NA, 1, 1), lwd = c(NA, 2, 2),
#        col = c("black", "blue"), pch = c(16, NA))
summary(Fit)
