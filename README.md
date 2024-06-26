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
