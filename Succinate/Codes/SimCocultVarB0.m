close all;
clear;

%This script is used to simulate co-culture experiment starting with
%different ratios or abundances. The last part tries to find a function
%describing the final observed proportions relative to initial
%conditions.
%Computes the line corresponding to the variety of attraction for different
%initial ratio. Uses the model without cross-feeding.

%Save or not the figure and data
save_data = 0; %If save = 1, figures and tables are saved (.fig or .xlsx) and overwrite the actual data. Otherwise, there are not.

%Load data
Name_Sheet = 'PVEPPU CO'; %Name of the sheet
data = readtable('Data/Ppu_Pve_growth_TECAN_succinate copy 2.xlsx', 'Sheet', Name_Sheet);
size_table = size(data);
Row_Start = 1;
Row_Fin = size_table(1);
%Time_step = table2array(data(Row_Start:Row_Fin, 2));
Time_step = 0:0.25:30;
data_Evol = table2array(data(Row_Start:Row_Fin, 3:(size_table(2)-2)));
CopyNum = 'VarietyConv'; 
%Pve values
load(strcat('./Data/','PVECDCWVal.mat')); 
mean_y_0_Pve = mean_y_0;
std_y_0_Pve = std_y_0;
table_Pve = load(strcat('./Data/','PVEKSLNTable.mat'));%Model @fun_Monod
mean_LN_mu_max_Pve = table_Pve.LN_mu_max(1);
var_LN_mu_max_Pve = table_Pve.LN_mu_max(2);
mean_LN_yield_Pve = table_Pve.LN_yield(1);
var_LN_yield_Pve = table_Pve.LN_yield(2);
mean_LN_Ks_Pve = table_Pve.LN_gamma(1);
var_LN_Ks_Pve = table_Pve.LN_gamma(2);
%Ppu values
load(strcat('./Data/','PPUCDCWVal.mat'));
mean_y_0_Ppu = mean_y_0;
std_y_0_Ppu = std_y_0;
table_Ppu = load(strcat('./Data/','PPUKSLNTable.mat'));
mean_LN_mu_max_Ppu = table_Ppu.LN_mu_max(1);
var_LN_mu_max_Ppu = table_Ppu.LN_mu_max(2);
mean_LN_yield_Ppu = table_Ppu.LN_yield(1);
var_LN_yield_Ppu = table_Ppu.LN_yield(2);
mean_LN_Ks_Ppu = table_Ppu.LN_gamma(1);
var_LN_Ks_Ppu = table_Ppu.LN_gamma(2);
%Parameters for the distribution generating the sample of X_m. Values drawn
%from Excel file sheet "Flow cytrometry"
ResInit = 5; %Initial resource concentration in mM.
mean_C_0 = 2.4*10^(-4);%ResInit*10^(-6)*48;
sigma_X_m = 0;
n_Time_step = length(Time_step);
props_pve = [0.1; 0.2; 0.5; 0.7; 0.8; 0.9; 0.995; 0.999];
mean_y_0_tot = mean_y_0_Pve;

%Output parameters
Num_Rep = 1;
tspan = [0 max(Time_step)]; %Time in hours
n_iter = 1;%5000;
Props = zeros(Num_Rep, length(props_pve));
num_fig = 1;

%The models
fun_Logistic = @(t, z, X_m, mu_max, yield, Ks) [z(1)*mu_max(1)*(1 - ((1/yield(1))*z(1) + (1/yield(2))*z(2))/X_m); z(2)*mu_max(2)*(1 - ((1/yield(1))*z(1) + (1/yield(2))*z(2))/X_m)];
fun_Monod = @(t, z, X_m, mu_max, yield, rho) [z(1)*mu_max(1)*((X_m - 1/yield(1)*z(1) - 1/yield(2)*z(2))/(X_m - 1/yield(1)*z(1) - 1/yield(2)*z(2) + rho(1)*X_m)); z(2)*mu_max(2)*((X_m - 1/yield(1)*z(1) - 1/yield(2)*z(2))/(X_m - 1/yield(1)*z(1) - 1/yield(2)*z(2) + rho(2)*X_m))];
% fun_Monod = @(t, z, X_m, mu_max, yield, Ks) [z(1)*mu_max(1)*((X_m - 1/yield(1)*z(1) - 1/yield(2)*z(2))/(X_m - 1/yield(1)*z(1) - 1/yield(2)*z(2) + Ks(1)*(2.4*10^(-4) + 1/yield(1)*2.96*10^(-7))));...
%     z(2)*mu_max(2)*((X_m - 1/yield(1)*z(1) - 1/yield(2)*z(2))/(X_m - 1/yield(1)*z(1) - 1/yield(2)*z(2) + Ks(2)*(2.4*10^(-4) + 1/yield(2)*2.78*10^(-7))))];

opts_1 = odeset('RelTol',1e-12,'AbsTol',1e-13,'NonNegative',1:2); %To smooth the curves obtained using ode45.

%Parameters generation. We want the same initial conditions for each ratio.
Mu_max_PVE =  lognrnd(mean_LN_mu_max_Pve, sqrt(var_LN_mu_max_Pve), 1, n_iter);
Mu_max_PPU =  lognrnd(mean_LN_mu_max_Ppu, sqrt(var_LN_mu_max_Ppu), 1, n_iter);
X_m = normrnd(mean_C_0, sigma_X_m, 1, n_iter);
yield_PVE =  lognrnd(mean_LN_yield_Pve, sqrt(var_LN_yield_Pve), 1, n_iter);
yield_PPU =  lognrnd(mean_LN_yield_Ppu, sqrt(var_LN_yield_Ppu), 1, n_iter);
Ks_PVE = lognrnd(mean_LN_Ks_Pve, 0*sqrt(var_LN_Ks_Pve), 1, n_iter);
Ks_PPU = lognrnd(mean_LN_Ks_Ppu, 0*sqrt(var_LN_Ks_Ppu), 1, n_iter);
yield = [yield_PVE yield_PPU];
mu_max = [Mu_max_PVE Mu_max_PPU];
Ks = [Ks_PVE Ks_PPU];

PVE_sim = zeros(length(Time_step), length(props_pve));
PPU_sim = zeros(length(Time_step), length(props_pve));
X_m_temp_tot = [];
for i = 1:length(props_pve)
    y_0_PVE = props_pve(i)*normrnd(mean_y_0_tot, 0*std_y_0_Pve, 1, n_iter);
    y_0_PPU = (1 - props_pve(i))*normrnd(mean_y_0_tot, 0*std_y_0_Ppu, 1, n_iter);
    Tab_output = {};
    for j = 1:n_iter
        mat_y_0 = [y_0_PVE(j) y_0_PPU(j)];
        %yield = yield + normrnd(0,0.1,1,2);
        X_m_temp = X_m(j) + (1/yield(1))*mat_y_0(1) + (1/yield(2))*mat_y_0(2);
        X_m_temp_tot(i) = X_m_temp;
        %System with 2 species
%         mu_max = mu_max + normrnd(0,0.1,1,2);
%         yield_test = [yield(1)+normrnd(0,0.1,1,2) yield(2)];
%         sol = ode45(@(t, y) fun_Logistic(t, y, X_m_temp, mu_max, yield, Ks), tspan,  mat_y_0, opts_1);
        sol = ode45(@(t, y) fun_Monod(t, y, X_m_temp, mu_max, yield, Ks), tspan,  mat_y_0, opts_1);
        z = deval(sol, Time_step);
        Tab_output{j} = [Time_step; z(1:2, :)];
    end
    %Randomly select nb_sim values among the n_iter values simulated
    ind_sim = randperm(n_iter);
    ind_sim = ind_sim(1:Num_Rep);

    for j = 1:Num_Rep
        Temp = Tab_output{ind_sim(j)};
        Temp = Temp';
        Props(j,i) = Temp(n_Time_step,2)/(Temp(n_Time_step,2) + Temp(n_Time_step,3));
        PVE_sim(:,i) = Temp(:,2);
        PPU_sim(:,i) = Temp(:,3);
    end
end

num_fig = num_fig + 1;
figure(num_fig)
plot(PVE_sim, PPU_sim, '--o')
hold on
fplot(@(x) yield(2)*X_m_temp - yield(2)/yield(1)*x, 'k-', [-5*10^(-5) 5*10^(-5)])
xlim([0 5*10^(-5)])
xlabel('Pve abundance')
ylabel('Ppu abundance')
legend(string(props_pve))
title('Evolution abundances')
num_fig = num_fig + 1;

% hold on
% plot(PVE_sim(1,:), PPU_sim(1,:), '-o')
% hold on
% fplot(@(x) yield(2)*(min(X_m_temp_tot) - X_m) - yield(2)/yield(1)*x, 'k-', [-5*10^(-5) 5*10^(-5)])

x_1 = SteadyState(PPU_sim(end, end), mat_y_0, mu_max);

if save_data == 1
    FolderName = strcat(cd, '/figures/');
    FigList = findobj(allchild(0), 'flat', 'Type', 'figure');
    for iFig = 1:length(FigList)
        FigHandle = FigList(iFig);
        FigName   = num2str(get(FigHandle, 'Number'));
        FigName = strcat('Simulated data coculture', FigName, CopyNum, num2str(ResInit), 'mM');
        set(0, 'CurrentFigure', FigHandle);
        savefig(fullfile(FolderName, [FigName '.fig']));
    end
end

%% Functions


function x_1 = SteadyState(x_2, x_0, mu_max)
x_1 = x_0(1)/(x_0(2)^(mu_max(1)/mu_max(2)))*x_2^(mu_max(1)/mu_max(2));
end