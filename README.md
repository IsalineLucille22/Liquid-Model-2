# Simulation of bacterial growth in a liquid suspension

## Overview

The different scripts are created in order to simulate mono or co-culture of Pseudomonas veronii and Pseudomonas putida in an homogeneous liquid suspension. With these scripts, you can especially test the effect of initial abundances, growth kinetics parameters and cross-feeding on the stationary state. 
The repository contains two folders. One folder is for experiments with a single resource, succinate, which is consumable by both species, leading to competition. The second folder contains data and scripts for experiments with two resources, putrescine and D-mannitol, each of which can be consumed by only one species.

## Link to BioRXiv

https://www.biorxiv.org/content/10.1101/2023.02.09.527847v1


## Succinate

This folder contains MATLAB codes for simulating our model, and a "Fitting" folder with MATLAB and R scripts used to fit our model parameters. Growth kinetics associted to each species are fitted on monocultures using the script named "ParamInfv2.m". Parameters corresponding to cross-feeding on byproducts are fitted on co-cultures using the R script named "ThresholdFitting.R". The last script in the folder is designed to compare simulated data to observations.
