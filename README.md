# Simulation of bacterial growth in a liquid suspension

## Overview

The different scripts are created in order to simulate mono or co-culture of Pseudomonas veronii and Pseudomonas putida in an homogeneous liquid suspension. With these scripts, you can especially test the effect of initial abundances, growth kinetics parameters and cross-feeding on the stationary state. 
The repository contains two folders. One folder is for experiments with a single resource, succinate, which is consumable by both species, leading to competition. The second folder contains data and scripts for experiments with two resources, putrescine and D-mannitol, each of which can be consumed by only one species.

## Link to BioRXiv

https://www.biorxiv.org/content/10.1101/2023.02.09.527847v1


## Succinate

This folder contains MATLAB codes for simulating our model, and a "Fitting" folder with MATLAB and R scripts used to fit our model parameters. The CoCultSimv2, SimCocultVarB0 and Main_Script are the script to be run to simulate the experiments. Parameters initialization for the different model are defined in these scripts.

**SimCocultVarB0**: Simulations of co-culture experiments starting with different initial ratios or abundances for each species using a model without cross-feeding. The only interactions are due to competition for the resource based on individual growth kinetics. Part of the script compares steady-state predictions with simulated time series produced by the model.

**CoCultSimv2**: This script simulates co-culture or monoculture experiments with the entire time-step evolution measured (10 replicates and 96 time steps). It is designed to compare simulations with observed data that include monocultures of PVE and PPU, as well as co-cultures with only a 1:1 ratio between the initial biomass of the two species. To simulate and compare other initial ratios between the two species, use the script named "Main_Script".

**Fitting folder**: Growth kinetics associted to each species are fitted on monocultures using the script named "ParamInfv2.m". Parameters corresponding to cross-feeding on byproducts are fitted on co-cultures using the R script named "ThresholdFitting.R". The last script in the folder is designed to compare simulated biomass to observed biomass.

## Models description

**Coculture** In co-cultures, in addition to the consumption of the primary resource according to intrinsic growth kinetics, we assume the possibility of cross-feeding or inhibition. Several functions can be chosen according to the desired model to simulate. They are,
1) fun_Hill_HandlingTime, Model: S_i + R -> P_i -> 2S_i, P_i <-> S_i + W_i, S_j + W_i -> P_j. Model with different wastes produced by each species. The waste produced by the other species can be used when its concentration is above a certain threshold. The threshold depends on the species and the detection function takes indicative form: max(z(j) - threshold(i), 0)/(z(j) - threshold(i))
2) fun_Inhibition, Model: S_i + R -> P_i -> 2S_i, P_i <-> S_i + W_i, S_j + W_i -> F. Model with different wastes produced by each species. The waste produced by the other species inhibits the development of the species consuming it when its concentration is above a certain threshold. The threshold depends on the species and the detection function takes indicative form: max(z(j) - threshold(i), 0)/(z(j) - threshold(i))
3) fun_Monod_tot_Death_rate, Model: S_i + R -> P_i -> 2S_i, P_i <-> S_i + W_i, S_j + W_i -> P_j, S_i -> 0. The last reaction corresponds to cell death rate (alpha). Model with different wastes produced by each species. The waste produced by the other species can be used when its concentration is above a certain threshold. The threshold depends on the species and the detection function takes indicative form: max(z(j) - threshold(i), 0)/(z(j) - threshold(i))
4) fun_Monod_tot, Model: S_i + R -> P_i -> 2S_i, P_i <-> S_i + W_i, S_j + W_i -> P_j. Model with different wastes produced by each species. The waste produced by the other species can be used when its concentration is above a certain threshold. The threshold depends on the species and the detection function takes indicative form: max(z(j) - threshold(i), 0)/(z(j) - threshold(i))
5) fun_Type_III, Model: S_i + R -> P_i -> 2S_i, P_i <-> S_i + W_i, S_j + W_i -> P_j.
6) fun_Unique_Waste, Model: S_i + R -> P_i -> 2S_i, P_i <-> S_i + W, S_j + W -> P_j. Model with one unique waste produced by each species. The waste produced can be used by both species with a different rate, when its concentration is above a certain threshold. The threshold does not depend on the species and the detection function takes indicative form: max(z(j) - threshold(i), 0)/(z(j) - threshold(i))

The biomass are defined as,
1) z(1), dx_1: biomass bacterial species 1.
2) z(2), dx_2: biomass bacterial species 2.
3) z(3), dy_1: biomass complex.
4) z(4), dy_2: biomass complex.
5) z(5), dw_1 biomass waste produced by species 1 that can be used by species 2.
6) z(6), dw_2: biomass waste produced by species 2 that can be used by species 1.
7) z(7), dr: biomass resource.

The choice of the model as well as the growth kinteics rates and interactions can be modified in the the "Main_Script.m" script.

**Monoculture** In a mono-culture with one nutrient, we assume species grow according to their growth kinetics by consuming the resource. A part of the consumed resource will be used for growth, while another part will be released as byproducts. The corresponding chemical reactions are,
S_i + R -> P_i -> 2S_i, P_i -> S_i + W_i

To simulate monoculture, the same functions as in co-cultures can be used; only the initial concentration of one of the two species should be set to 0.


