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

#Change it according to the data folder
setwd("/Users/isalinelucille-guex/Documents/Liquid-models/Succinate/Codes/Data")

Data_1 = read_excel("Ppu_Pve_toluene_succ_replicate_growth.xlsx", sheet = 3)

Ratio <- Data_1$Ratio
Data_Tania <- data.frame(
  Ratio = Ratio,
  x = c(Data_1$Abundance1, Data_1$Abundance2, Data_1$Abundance3, Data_1$Abundance4),
  y = c(Data_1$Abundance5, Data_1$Abundance6, Data_1$Abundance7, Data_1$Abundance8),
  z = c(Data_1$Abundance9, Data_1$Abundance10, Data_1$Abundance11, Data_1$Abundance12),
  v = c(Data_1$Abundance13, Data_1$Abundance14, Data_1$Abundance15, Data_1$Abundance16)
)
R_0 = 2.4*10^(-4)
TT = seq(0,24,0.5)
nb_ratio = length(Data_Tania$Ratio)

## =======================================================================
## Model double wastes fitting
## =======================================================================

fun_Monod_tot <- function(t, state, parms, parms_est, index){
  with(as.list(c(state, parms, parms_est)), {
    dx_1 = (2*kappa_12 + kappa_13)*y_1 - kappa_11*x_1*r - (w_2 - threshold_1 > 0)*kappa_14*x_1*w_2
    dx_2 = (2*kappa_22 + kappa_23)*y_2 - kappa_21*x_2*r - (w_1 - threshold_2 > 0)*kappa_24*x_2*w_1
    dy_1 = -(kappa_12 + kappa_13)*y_1 + kappa_11*x_1*r + (w_2 - threshold_1 > 0)*kappa_14*x_1*w_2
    dy_2 = -(kappa_22 + kappa_23)*y_2 + kappa_21*x_2*r + (w_1 - threshold_2 > 0)*kappa_24*x_2*w_1
    dw_1 = kappa_13*y_1 - (w_1 - threshold_2 > 0)*kappa_24*x_2*w_1
    dw_2 = kappa_23*y_2 - (w_2 - threshold_1 > 0)*kappa_14*x_1*w_2
    dr = -kappa_11*x_1*r - kappa_21*x_2*r;
    list(c(dx_1, dx_2, dy_1, dy_2, dw_1, dw_2, dr))
  }) 
}

fun_Monod_Hill <- function(t, state, parms, parms_est, index){
  with(as.list(c(state, parms, parms_est)), {
    Hill_Coeff_1 = 1/(1 + (threshold_1/w_2)^k_1);  Hill_Coeff_2 = 1/(1 + (threshold_2/w_1)^k_1);
    dx_1 = (2*kappa_12 + kappa_13)*y_1 - kappa_11*x_1*r - Hill_Coeff_1*kappa_14*x_1*w_2
    dx_2 = (2*kappa_22 + kappa_23)*y_2 - kappa_21*x_2*r - Hill_Coeff_2*kappa_24*x_2*w_1
    dy_1 = -(kappa_12 + kappa_13)*y_1 + kappa_11*x_1*r + Hill_Coeff_1*kappa_14*x_1*w_2
    dy_2 = -(kappa_22 + kappa_23)*y_2 + kappa_21*x_2*r + Hill_Coeff_2*kappa_24*x_2*w_1
    dw_1 = kappa_13*y_1 - Hill_Coeff_2*kappa_24*x_2*w_1
    dw_2 = kappa_23*y_2 - Hill_Coeff_1*kappa_14*x_1*w_2
    dr = -kappa_11*x_1*r - kappa_21*x_2*r;
    list(c(dx_1, dx_2, dy_1, dy_2, dw_1, dw_2, dr))
  }) 
}

fun_Type_III <- function(t, state, parms, parms_est, index){
  with(as.list(c(state, parms, parms_est)), {
    Hill_Coeff_1 = k_1;  Hill_Coeff_2 = k_2;
    dx_1 = (2*kappa_12 + kappa_13)*y_1 - kappa_11*x_1*r - kappa_14*x_1*w_2^Hill_Coeff_1
    dx_2 = (2*kappa_22 + kappa_23)*y_2 - kappa_21*x_2*r - kappa_24*x_2*w_1^Hill_Coeff_2
    dy_1 = -(kappa_12 + kappa_13)*y_1 + kappa_11*x_1*r + kappa_14*x_1*w_2^Hill_Coeff_1
    dy_2 = -(kappa_22 + kappa_23)*y_2 + kappa_21*x_2*r + kappa_24*x_2*w_1^Hill_Coeff_2
    dw_1 = kappa_13*y_1 - kappa_24*x_2*w_1^Hill_Coeff_2
    dw_2 = kappa_23*y_2 - kappa_14*x_1*w_2^Hill_Coeff_1
    dr = -kappa_11*x_1*r - kappa_21*x_2*r;
    list(c(dx_1, dx_2, dy_1, dy_2, dw_1, dw_2, dr))
  }) 
}

ModelCost <- function(P) {
  diff = numeric()
  for(i in 1:nb_ratio){
    x_10 = Data_Tania[i,3] #x_0 Pve
    x_20 = Data_Tania[i,4] #x_0 Ppu
    parameters_temp = c(parameters, x_10, x_20)
    state = c(x_1 = x_10, x_2 = x_20, y_1 = 0, y_2 = 0, w_1 = 0, w_2 = 0, r = 2.4*10^(-4))
    out <- ode(y = state, func = fun_Monod_Hill, parms = parameters_temp,  parms_est = P, index = i, times = TT, atol = 1e-9, rtol = 1e-9)
    model = out[length(out[,1]),]
    Ppu_fin = model[3] + model[5] #Final Ppu biomass
    model = model[2] + model[4] #Final Pve biomass
    # diff = c(diff, abs(Data_Tania[i,2] - model)) #For Absolute biomass
    # diff = c(diff, abs(Data_Tania[i,2] - model) + abs(Data_Tania[i,5] - Ppu_fin)) #For Absolute biomass considering Pve and Ppu
    diff = c(diff, abs((Data_Tania[i,2]/(Data_Tania[i,2] + Data_Tania[i,5])) - model/(model + Ppu_fin)) + abs(Data_Tania[i,2] - model) + abs(Data_Tania[i,5] - Ppu_fin)) #For ratio
  }
  Data_Tania_2 = Data_Tania[,1:2]
  Data_Tania_2[,2] = diff
  return(diff)
}

parameters <- c(kappa_11 = 250900, kappa_12 = 0.48, kappa_13 = 1.26, kappa_14 = 250900, kappa_21 = 200000, kappa_22 = 1.09, kappa_23 = 1.65, kappa_24 = 200000)
parameters_init <- c(threshold_1 = 2.7e-04, threshold_2 = 2.7e-04, k_1 = 5)
# parameters_init <- c(threshold_1 = 2.7e-04, threshold_2 = 2.7e-04, k_1 = 5, k_2 = 5)



x_10 = Data_Tania[1,3]
x_20 = Data_Tania[1,4]
Fit <- modFit(f = ModelCost, p = parameters_init, lower = rep(0, 3),
              upper = c(1e-03, 1e-03, 50))

# Fit <- modFit(f = ModelCost, p = parameters_init, lower = rep(0, 4),
#               upper = c(1e-03, 1e-03, 50, 50))




