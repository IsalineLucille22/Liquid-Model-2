clear;
close all;

%Model considered are of type: S + R -> P -> 2S, P <-> S + W, P -> S + F.
%Model with double waste, one (W) that can be reused by both species and the
%other (F) can't be reused.

% Enter parameters to test:

Name_file = 'test'; %Name the produced data.
save_data = 0; %If save = 1, figures and tables are saved (.fig/.pdf or .xlsx) and overwrite the actual file of the same name.

%% Allocation of the parameter values

%Parameters allocated accoridng to observed values
table_name = ['PVEManuAdaptScalLN'; 'PPUManuAdaptScalLN'];
%Estimated values
%Pve values
table_Pve = load(strcat('./Data/',table_name(1,:), '.mat'));
mean_LN_k1_Pve = table_Pve.LN_k1; %Kappa_1
mean_LN_k2_Pve = table_Pve.LN_k2; %Mu_max
mean_LN_k3_Pve = table_Pve.LN_k3; %Golbal yield
%Ppu values
table_Ppu = load(strcat('./Data/',table_name(2,:), '.mat'));
mean_LN_k1_Ppu = table_Ppu.LN_k1; %Kappa_1
mean_LN_k2_Ppu = table_Ppu.LN_k2; %Mu_max
mean_LN_k3_Ppu = table_Ppu.LN_k3; %Global yield

%If set mu_max to empty vector, mu_max = [], then estimated values from
%parameters inference are used (Table P__ManuAdptScalLN.mat) to generate kappa_2 (mu_max) and kappa_3
%(global yield). Values have to be given to the other parameters.
mu_max = [];%[0.48 1.05];%[0.52 1.2]; %Maximum growth rates of both species
std_mu_max = [0 0]; %Standard deviation mu_max. 

%The rate kappa_1 (h*g/ml)^(-1). They should be
%of the order 1/K_s where K_s is the half velocity constant in mL.
kappa_1 = [2.5090e+05; 2e+05]; %Rates for reaction S_i + R -> P_i
% where S_i is the species, R the unique resource initially present and P_i the complex corresponding to species i.
std_kappa_1 = [0 0]; %Standard deviation kappa_1.

%Concentration threshold. When the concentration of the reusable waste (W)
%is above this threshold, then the species i can use it.
% Threshold = 0*[1.400e-04; 2.8000e-05]; %Threshold values
%100:1, 10:1, 1:1, 1:10, 1:100
%T(1,:) for used of W_2 by PVE, T(2,:) for used of W_1 by PPU
%Indicative function
% Threshold = [1.397e-04 1.397e-04 1.397e-04 1.397e-04 1.4345e-04; 6.5e-05 1.0e-05 2.8e-05 2.8e-05 2.8e-05];
% Threshold = [1.43e-04 1.43e-04 1.43e-04 1.43e-04 1.43e-04; 2.8e-05 2.8e-05 2.8e-05 2.8e-05 2.8e-05];
%Sigmoid function
Threshold = [1.726828e-04 1.726828e-04 1.726828e-04 1.726828e-04 1.726828e-04; 4.728681e-05 4.728681e-05 4.728681e-05 4.728681e-05 4.728681e-05];%
% Inhibition function
% Threshold = [8.854093e-05 8.854093e-05 8.854093e-05 8.854093e-05 8.854093e-05; 2.189977e-04 2.189977e-04 2.189977e-04 .189977e-04 2.189977e-04];% Bilateral inhibition
% Threshold = [2.700000e-04 2.700000e-04 2.700000e-04 2.700000e-04 2.700000e-04; 9.029489e-04 9.029489e-04 9.029489e-04 9.029489e-04 9.029489e-04];% Unilateral inhibition

mean_R_0 = 2.4*10^(-4); %Initial resource concentration
sigma_R_0 = 0*0.5*10^(-5); %Standard deviation for the resource concentration.

Time_step = 0:0.25:24; %Time step, 0:0.25:24 corresponds to the real measurement times of the experiment.
%% Test 5 ratios
close all

Name_Model = @fun_Hill_HandlingTime; %Type of the function to test. Description in the corresponding function file.

num_ratio_min = 1; %Number of the ratio tested. In experiment, there are 5 different ratios 100:1, 10:1, 1:1, 1:10, 1:100.
num_ratio_max = 5;
ratios = [num_ratio_min num_ratio_max];

[Fin_Abund_PVE, Fin_Abund_PPU, kappa_mat, mean_PVE_fin_props, Props_sim_PVE, Props_obs_PVE] = Main_Fun(Name_Model, mu_max, std_mu_max, kappa_1, std_kappa_1,...
    mean_LN_k1_Pve, mean_LN_k2_Pve, mean_LN_k3_Pve, mean_LN_k1_Ppu, mean_LN_k2_Ppu, mean_LN_k3_Ppu,...
    Threshold, ratios, mean_R_0, sigma_R_0, Time_step, Name_file, save_data);