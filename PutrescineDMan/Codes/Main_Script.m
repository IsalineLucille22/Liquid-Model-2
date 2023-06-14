clear;
close all;

%Model considered are of type: S_i + R_i -> P_i -> 2S_i, P_i <-> S_i + W_i
%Model without cross-feeding is generated using function
%fun_Hill_HandlingTime. Model with cross-feeding is generated using
%function fun_Hill_HandlingTimev3. Threshold value can be changed
%into Main_Fun.m

% Enter parameters to test:

Name_file = 'FigureCo'; %Name the produced data.
save_data = 0; %If save = 1, figures and tables are saved (.fig/.pdf or .xlsx) and overwrite the actual file of the same name.

%Load observed data from Excel files
Name_Sheet = [{'CoCult (PVE PPU 100_1)'}, {'CoCult (PVE PPU 10_1)'}, {'CoCult (PVE PPU 1_1)'},{'CoCult (PVE PPU 1_10)'}, {'CoCult (PVE PPU 1_100)'}]; 
%For Fluo [{'CoCult (PVE PPU 100_1 Fluo)'}, {'CoCult (PVE PPU 10_1 Fluo)'}, {'CoCult (PVE PPU 1_1 Fluo)'},{'CoCult (PVE PPU 1_10 Fluo)'}, {'CoCult (PVE PPU 1_100 Fluo)'}];%
num_ratio = 1:5;
DataFile = strcat('./Data/', Name_file, '.mat'); 

%Estimated values
%Pve values
table_Pve = load(strcat('./Data/', 'PVED-ManLNTablev2.mat'));
mu_max_Pve = table_Pve.LN_k2; %Mu_max
mean_LN_k1_Pve = table_Pve.LN_k1; %kappa_1
mean_LN_k3_Pve = table_Pve.LN_k3; %Golbal yield
table_Pve = load(strcat('./Data/', 'PVED-ManPutresLNTable.mat'));
mean_LN_k2_Pve = table_Pve.LN_k2; %Mu_max
%Ppu values
table_Ppu = load(strcat('./Data/', 'PPUPutresLNTablev2.mat'));
mean_LN_k1_Ppu = table_Ppu.LN_k1; %kappa_1
mean_LN_k2_Ppu = table_Ppu.LN_k2; %Mu_max
%     table_Ppu = load(strcat('./Data/', 'PPUD-ManPutresLNTable.mat'));
%     mu_max_Ppu = table_Ppu.LN_k2; %Mu_max
mean_LN_k3_Ppu = table_Ppu.LN_k3; %Global yield 
yield_Ppu = 0.34; %Based on histogram 3PPUD-ManPutres. Logistic estimation

%     100:1, 10:1, 1:1, 1:10, 1:100
%     Threshold values
Threshold_values = [3.4e-04 2.85e-04 2.4e-04 1.8e-04 1.2e-04;... %Use W_2 by S_1
    12.9e-04 12.9e-04 12.9e-04 12.9e-04 12.9e-04;... %Use W_1 by S_2
    7.35e-04 12.5e-04 6.35e-04 15.8e-04 14e-04;]; %Use W_1 by S_1

Mat_y_0 = [4.94e-06 2.29e-08; 4.9e-06 1.75e-07; 2.44e-06 1.01e-06;...
    5.16e-07 1.84e-06; 7.36e-08 1.78e-06];%Initial abundances measured. First row for P.veronni. Second row for P.putida.

mean_R_0 = [6.67*10^(-3)*72/1000; 1*10*10^(-3)*48/1000]; %Initial resource concentration, correspond to 10mM of succinate. First component for D-Mannitol, second component for Putrescine.
sigma_R_0 = [0; 0]; %Standard deviation for the resource concentration. When 0, there is no variation in the resource concentrations. First component for D-Mannitol, second component for Putrescine.

Time_step = 0:0.5:48; %Time step, it corresponds to the real measurement times of the experiment.

Name_Model = @fun_Hill_HandlingTimev3; %Type of the function to test. Description in the corresponding function file.

[Tot_Gram_Evol, Abund_sim_PVE, Abund_sim_PPU, kappa_mat] = Main_Fun(Name_Model, mean_R_0, sigma_R_0, Mat_y_0, Time_step, Name_file, Name_Sheet, num_ratio, DataFile, mu_max_Pve, mean_LN_k2_Pve, mean_LN_k1_Pve, mean_LN_k3_Pve, mean_LN_k1_Ppu, mean_LN_k2_Ppu, mean_LN_k3_Ppu, yield_Ppu, Threshold_values, save_data);