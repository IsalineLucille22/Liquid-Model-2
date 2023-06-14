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
# setwd("/Users/isalinelucille-guex/switchdrive/shared UNIL/Ppu_Pve_growth_data/version 2/Maxime data/Dmannitol or Putrescine or No Carbon")
setwd("/Users/guexi/switchdrive/shared UNIL/Ppu_Pve_growth_data/version 2/Maxime data/Dmannitol or Putrescine or No Carbon")

Data_1 = read_excel("DataMaximev4.xlsx", sheet = 3)

Ratio <- Data_1$`Ratio (Pve:Ppu)`
Data_Tania <- data.frame(
  Ratio = Ratio,
  x = c(Data_1$Abundance1),
  y = c(Data_1$Abundance5),
  z = c(Data_1$Abundance9),
  v = c(Data_1$Abundance13)
)
R_0 = c(6.67*10^(-3)*72/1000, 1*10*10^(-3)*48/1000)
TT = seq(0,48,0.5)
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
    Hill_Coeff_1 = 1/(1 + (threshold_1/w_2)^k_1);  Hill_Coeff_2 = 0*1/(1 + (threshold_2/w_1)^k_2);
    dx_1 = (2*kappa_12 + kappa_13)*y_1 - kappa_11*x_1*r_1 - Hill_Coeff_1*kappa_14*x_1*w_2
    dx_2 = (2*kappa_22 + kappa_23)*y_2 - kappa_21*x_2*r_2 - Hill_Coeff_2*kappa_24*x_2*w_1
    dy_1 = -(kappa_12 + kappa_13)*y_1 + kappa_11*x_1*r_1 + Hill_Coeff_1*kappa_14*x_1*w_2
    dy_2 = -(kappa_22 + kappa_23)*y_2 + kappa_21*x_2*r_2 + Hill_Coeff_2*kappa_24*x_2*w_1
    dw_1 = kappa_13*y_1 - Hill_Coeff_2*kappa_24*x_2*w_1
    dw_2 = kappa_23*y_2 - Hill_Coeff_1*kappa_14*x_1*w_2
    dr_1 = -kappa_11*x_1*r_1 
    dr_2 = - kappa_21*x_2*r_2;
    list(c(dx_1, dx_2, dy_1, dy_2, dw_1, dw_2, dr_1, dr_2))
  }) 
}

fun_Monod_Hill_2 <- function(t, state, parms, parms_est, index){
  with(as.list(c(state, parms, parms_est)), {
    # Hill_Coeff_1 = 1/(1 + (threshold_1/(x_1/x_2))^k_1);  Hill_Coeff_2 = 1/(1 + (threshold_2/(x_1/x_2))^k_2);
    Hill_Coeff_1 = (1 - 1/(1 + (threshold_1/parms[9])^k_1));  Hill_Coeff_2 = 0*1/(1 + (threshold_2/parms[10])^k_2);
    dx_1 = (2*kappa_12 + kappa_13)*y_1 - kappa_11*x_1*r_1 - Hill_Coeff_1*kappa_14*x_1*w_2
    dx_2 = (2*kappa_22 + kappa_23)*y_2 - kappa_21*x_2*r_2 - Hill_Coeff_2*kappa_24*x_2*w_1
    dy_1 = -(kappa_12 + kappa_13)*y_1 + kappa_11*x_1*r_1 + Hill_Coeff_1*kappa_14*x_1*w_2
    dy_2 = -(kappa_22 + kappa_23)*y_2 + kappa_21*x_2*r_2 + Hill_Coeff_2*kappa_24*x_2*w_1
    dw_1 = kappa_13*y_1 - Hill_Coeff_2*kappa_24*x_2*w_1
    dw_2 = kappa_23*y_2 - Hill_Coeff_1*kappa_14*x_1*w_2
    dr_1 = -kappa_11*x_1*r_1 
    dr_2 = -kappa_21*x_2*r_2;
    list(c(dx_1, dx_2, dy_1, dy_2, dw_1, dw_2, dr_1, dr_2))
  }) 
}

ModelCost <- function(P) {
  diff = numeric()
  for(i in 1:nb_ratio){
    x_10 = Data_Tania[i,3] #x_0 Pve
    x_20 = Data_Tania[i,4] #x_0 Ppu
    state = c(x_1 = x_10, x_2 = x_20, y_1 = 0, y_2 = 0, w_1 = 0, w_2 = 0, r_1 = 0.00048024, r_2 = 0.00048)
    parameters_temp = c(parameters, x_10, x_20)
    # print(parameters_temp[9])
    out <- ode(y = state, func = fun_Monod_Hill_2, parms = parameters_temp,  parms_est = P, index = i, times = TT, atol = 1e-9, rtol = 1e-9)
    model = out[length(out[,1]),]
    Ppu_fin = model[3] + model[5] #Final Ppu biomass
    model = model[2] + model[4] #Final Pve biomass
    diff = c(diff, Data_Tania[i,2] - model) #For Absolute biomass
    # diff = c(diff, abs(Data_Tania[i,2] - model) + abs(Data_Tania[i,5] - Ppu_fin)) #For Absolute biomass considering Pve and Ppu
    # diff = c(diff, abs((Data_Tania[i,2]/(Data_Tania[i,2] + Data_Tania[i,5])) - model/(model + Ppu_fin)) + abs(Data_Tania[i,2] - model))# + abs(Data_Tania[i,5] - Ppu_fin)) #For ratio
  }
  Data_Tania_2 = Data_Tania[,1:2]
  Data_Tania_2[,2] = diff
  return(diff)
}

parameters <- c(kappa_11 = 4.6589e+03, kappa_12 = 0.2718, kappa_13 = 0.8898, kappa_14 = 4.6589e+03, kappa_21 = 1.4235e+04, kappa_22 = 1.5386, kappa_23 = 4.2391, kappa_24 = 1.4235e+04)
parameters_init <- c(threshold_1 = 2.7e-08, threshold_2 = 2.7e-06, k_1 = 5, k_2 = 5)


x_10 = Data_Tania[1,3]
x_20 = Data_Tania[1,4]
Fit <- modFit(f = ModelCost, p = parameters_init, lower = rep(0, 4),
              upper = c(1e-03, 1e-03, 50, 150))
# Fit <- modFit(f = ModelCost, p = parameters_init, lower = rep(0, 4),
#               upper = c(4.6589e+04, 4.6589e+04, 150, 150))