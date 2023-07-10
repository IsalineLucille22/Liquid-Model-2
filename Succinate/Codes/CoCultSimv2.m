clear;
close all;

%This script simulates the co-culture or mono-culture experiments with the whole time step
%evolution (10 replicates and 96 time steps). It is made to compare with
%the results of the file "Ppu_Pve_growth_TECAN_succinate copy 2" that
%contains mono-cultures of PVE and PPU (also used for fitting) and
%co-culture but only with 1:1 ratio. To compare the other ratios, use the
%file Main_Script. 

%Save or not the figure and data
save_data = 0; %If save = 1, figures and tables are saved (.fig or .xlsx) and overwrite the actual data. Otherwise, there are not.

%Load data
CopyNum = 'FluoFinalFluoValues'; %Name of the saved data
Name_Sheet = 'Flow cytometry'; % Name of the sheet
data_prop = readtable('Data/Ppu_Pve_growth_TECAN_succinate copy 2.xlsx', 'Sheet', Name_Sheet,'Format','auto');
Props_obs_PVE = str2double(table2array(data_prop(35:37,4)));%For fluo: 43:52,6 For OD: 35:37,4
Props_obs_PPU = str2double(table2array(data_prop(35:37,5)));%For fluo: 43:52,7 For OD: 35:37,5
PVE_gain_obs = str2double(table2array(data_prop(35:37,6)));
PPU_gain_obs = str2double(table2array(data_prop(35:37,7)));
Name_Sheet = 'PVEPPU CO'; %Name of the sheet 'PVEPPU CO', 'PPU alone', 'PVE alone', 'PVE alone (OD)', 'PPU alone (OD)'
%Name_Sheet = ['PVEPPU CO PVE (Fluo)'; 'PVEPPU CO PPU (Fluo)'];
data_Evol = {};
for i = 1:length(Name_Sheet(:,1))
    data = readtable('Data/Ppu_Pve_growth_TECAN_succinate copy 2.xlsx', 'Sheet', Name_Sheet(i,:),'Format','auto');
    size_table = size(data);
    Row_Start = 1;
    Row_Fin = size_table(1);
    Time_step = table2array(data(Row_Start:Row_Fin, 2));
    data_Evol{i} = table2array(data(Row_Start:Row_Fin, 3:(size_table(2)-2)));
end


%Pve values
load('Data/PVECDCWVal.mat'); 
mean_y_0_Pve = mean_y_0;
std_y_0_Pve = std_y_0;
table_Pve = load('Data/PVEManuAdaptScalLN.mat');
mean_LN_k1_Pve = table_Pve.LN_k1; %Kappa_1
mean_LN_k2_Pve = table_Pve.LN_k2; %Mu_max
mean_LN_k3_Pve = table_Pve.LN_k3; %Golbal yield
%Ppu values
load('Data/PPUCDCWVal.mat');
mean_y_0_Ppu = mean_y_0;
std_y_0_Ppu = std_y_0;
table_Ppu = load('Data/PPUManuAdaptScalLN.mat');
% table_Ppu = load('PPULagLN.mat');
mean_LN_k1_Ppu = table_Ppu.LN_k1; %Kappa_1
mean_LN_k2_Ppu = table_Ppu.LN_k2; %Mu_max
mean_LN_k3_Ppu = table_Ppu.LN_k3; %Global yield


%Parameters for the distribution generating the sample of X_m. Values drawn
%from Excel file sheet "Flow cytrometry"
mean_R_0 = 2.4*10^(-4);
sigma_R_0 = 0;
n_Time_step = length(Time_step);

%Output parameters
Num_Rep = 3; %For fluo:10, For OD: 3
tspan = [0 max(Time_step)]; %Time in hours

opts_1 = odeset('RelTol',1e-8,'AbsTol',1e-9,'NonNegative',1:7); %To smooth the curves obtained using ode45.


%Parameters generation
n_iter = 500;
X_m = normrnd(mean_R_0, sigma_R_0, 1, n_iter); 
y_0_PVE = 1*1/2*normrnd(mean_y_0_Pve, std_y_0_Pve, 1, n_iter);
y_0_PPU = 1*1/2*normrnd(mean_y_0_Ppu, std_y_0_Ppu, 1, n_iter);
%Computation of the rates of the reactions
kappa_2 = [lognrnd(mean_LN_k2_Pve(1), mean_LN_k2_Pve(2),  1, n_iter);...
    lognrnd(mean_LN_k2_Ppu(1), mean_LN_k2_Ppu(2),  1, n_iter)];
kappa_3 = [lognrnd(mean_LN_k3_Pve(1), 0*mean_LN_k3_Pve(2),  1, n_iter);...
    lognrnd(mean_LN_k3_Ppu(1), 0*mean_LN_k3_Ppu(2),  1, n_iter)];
kappa_1 = [lognrnd(mean_LN_k1_Pve(1), 0*mean_LN_k1_Pve(2),  1, n_iter);...
    lognrnd(mean_LN_k1_Ppu(1), 0*mean_LN_k1_Ppu(2),  1, n_iter)];
% kappa_1 = [2.5090e+05 + normrnd(0, 0, 1, n_iter); 2.0e+05 + normrnd(0, 0, 1, n_iter)];

Tab_output = {};
for i = 1:n_iter
    mat_y_0 = [1*y_0_PVE(i) 1*y_0_PPU(i) 0 0 0 0 mean_R_0 0 0 0];
    kappa_mat = [kappa_1(1,i) kappa_2(1,i) kappa_3(1,i) kappa_1(1,i); kappa_1(2,i) kappa_2(2,i) kappa_3(2,i) kappa_1(2,i)];
    sol = ode45(@(t, y) fun_Monod_tot(t, y, kappa_mat, [1.397e-04; 2.8e-05]), tspan,  mat_y_0, opts_1);
%     sol = ode45(@(t, y) fun_Hill_HandlingTime(t, y, kappa_mat, [1.726828e-04; 4.728681e-05]), tspan,  mat_y_0, opts_1);
    z = deval(sol, Time_step);
    X_1 = z(1, :) + z(3, :);
    X_2 = z(2, :) + z(4, :);
    Tab_output{i} = [Time_step'; X_1; X_2];
end
%Randomly select nb_sim values among the n_iter values simulated
ind_sim = randperm(n_iter);
ind_sim = ind_sim(1:Num_Rep);


Tot_Fin_OD = zeros(1, Num_Rep);
Props = zeros(1, Num_Rep);
PVE_fin_abs = zeros(1, Num_Rep);
PPU_fin_abs = zeros(1, Num_Rep);
[Tot_OD_Evol,Pve_OD, Ppu_OD] = deal(zeros(length(Time_step), Num_Rep));
num_fig = 1;
for i = 1:Num_Rep
    Temp = Tab_output{ind_sim(i)};
    Temp = Temp';
    Pve_OD(:,i) = Temp(:,2)/(5.4318*10^(-4)); %For biomass in cell carbon dry weight
    Ppu_OD(:,i) = Temp(:,3)/(5.9484*10^(-4)); %For biomass in cell carbon dry weight
    PVE_fin_abs(i) = Temp(n_Time_step,2);
    PPU_fin_abs(i) = Temp(n_Time_step,3);
    Tot_Fin_OD(i) = Pve_OD(n_Time_step) + Ppu_OD(n_Time_step) + 0.086;
    Tot_OD_Evol(:,i) = Pve_OD(:,i) + Ppu_OD(:,i) + 0.086; %For OD
%     Tot_OD_Evol(:,i) = Pve_OD(:,i) + Ppu_OD(:,i); %For Fluo
    Props(i) = Temp(n_Time_step,2)/(Temp(n_Time_step,3)+Temp(n_Time_step,2));
end

[t_test, h_test] = ttest2(Props', Props_obs_PVE);

figure(num_fig)
time_end = length(Time_step);

Tot_OD_Evol = {Tot_OD_Evol}; %For OD
% Tot_OD_Evol = {Pve_OD + 0.086, Ppu_OD + 0.086}; %For Fluo
Time_Lag = [1; 1];%[25; 14];
for i = 1:length(Name_Sheet(:,1))
    std_Sim = std(Tot_OD_Evol{i},0,2);
    mean_Sim = mean(Tot_OD_Evol{i},2);
    mean_Test = mean(data_Evol{i},2);
    std_Test = std(data_Evol{i},0,2);
    errorbar(Time_step(1:time_end), mean_Test(1:time_end), std_Test(1:time_end),'r--o','MarkerSize',1,'LineWidth', 1,'MarkerEdgeColor','black','MarkerFaceColor','black');
%     errorbar(Time_step(1:time_end-(Time_Lag(i) - 1)), mean_Test(Time_Lag(i):time_end), std_Test(Time_Lag(i):time_end),'r--o','MarkerSize',1,'LineWidth', 1,'MarkerEdgeColor','black','MarkerFaceColor','black');
    hold on
    errorbar(Time_step(1:time_end), mean_Sim(1:time_end), std_Sim(1:time_end),'b--o','MarkerSize',1,'LineWidth', 1,'MarkerEdgeColor','red','MarkerFaceColor','red');
    hold on
end
xlabel('Time (hours)')
ylabel('Abundance (OD)')
title('Optical density evolution')
num_fig = num_fig + 1;

figure(num_fig)
boxplot([Props; Props_obs_PVE']', 'Labels',{'Sim','Obs'})
hold on 
scatter(ones(size(Props)).*(1+(rand(size(Props)) - 0.5)/10), Props, 'r', 'x')
scatter(ones(size(Props)).*(2+(rand(size(Props)) - 0.5)/10), Props_obs_PVE, 'r', 'x')
xlabel('Category')
ylabel('Proportions')
title('Boxplots final proportions PVE');

if save_data == 1
    FolderName = strcat(cd, '/figures/');
    FigList = findobj(allchild(0), 'flat', 'Type', 'figure');
    for iFig = 1:length(FigList)
        FigHandle = FigList(iFig);
        FigName = num2str(get(FigHandle, 'Number'));
        FigName = strcat(FolderName, 'Simudataco', FigName, CopyNum);
        set(0, 'CurrentFigure', FigHandle);
        saveas(FigHandle, FigName, 'pdf');
    end
end