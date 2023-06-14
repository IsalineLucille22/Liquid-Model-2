clear;
close all;

%This script is used to fit parameters for Maxime's data

%Save or not the figure and data
save_data = 0; %If save = 1, figures and tables are saved (.fig or .xlsx) and overwrite the actual data. Otherwise, there are not.

%Load data
Name_Sheet = 'Mono PVE Putres g per mL';%'Mono PPU (g per mL)';%'Mono PPU Putres g per mL';%'Mono PVE D-Man g per mL'; %
data = readtable('../Data/DataMaximev4.xlsx', 'Sheet', Name_Sheet);%readtable('../Data/DataMaximev3.xlsx', 'Sheet', Name_Sheet);
CopyNum = 'PVE-Putres';%'PVED-Man'; %'PPUPutres'
DataFile = strcat('./Data/', Name_Sheet, CopyNum, '.mat');
size_table = size(data);
Row_Start = 1;
Row_Fin = size_table(1);%40;%size_table(1); 
Time_step = table2array(data(Row_Start:Row_Fin, 1)) - table2array(data(Row_Start, 1)); %Time of the observations in hours
data_Evol = table2array(data(1:size_table(1), 4:(size_table(2))));
data_Evol(data_Evol<0) = 0;
mean_R_0 = 10*10^(-3)*48/1000;%+ 6.67*10^(-3)*72/1000;%
sigma_R_0 = 0;
Num_Rep = length(data_Evol(1,:));
data_Train = data_Evol(Row_Start:Row_Fin, :);%data_Evol(Row_Start:Row_Fin, (mod(1:Num_Rep,2) == 1));%(mod(1:Num_Rep,1) == 1)); %Select all data
data_Test = data_Evol(Row_Start:Row_Fin, :);%data_Evol(Row_Start:Row_Fin, (mod(1:Num_Rep,2) == 1)); %Select half the data to compare it with simulation (won't be used before section 4)
num_fig = 1;
gamma = std(data_Evol(Row_Start:Row_Fin,:),0,2);%2*10^(-5)*ones(Row_Fin - Row_Start + 1,1);%
gamma = gamma';

opt_1 = odeset('RelTol',1e-9,'AbsTol',1e-9);%,'NonNegative',1);
opt_2 = odeset('RelTol',1e-9,'AbsTol',1e-9);%, 'NonNegative',1:4); %To smooth the curves obtained using ode45.

%% Mu_max, yield inference
%Computation of the parameters used into the model (Logistic for the estimate of mu_max and yield). Using a MH algorithm,
%the goal of this section is to generate samples for parameters mu_max and
%X_m. These sample will be then used to determine log-normal priors for our model. There are generated using half the data (training set)
n_time = length(Time_step);

tspan = [min(Time_step), max(Time_step)]; %Time interval in hours

pos_guess = [2 3];
nb_param = length(pos_guess);
vect_alpha = [1; 0; 0; 0];%Assuming the parameters follow a uninfromative prior (uniform);
vect_beta = [1; 15; 10; 1];
T_sample = 5000; %Length of the samples
T_0_sample = 50; %Remove first values
T_fin_sample = T_sample - T_0_sample; %Observations kept in order to form samples
R_0_temp = normrnd(mean_R_0, sigma_R_0);
mean_param = zeros(1,nb_param);
vect_init = [R_0_temp; 0.01; 0.5; 0.1*R_0_temp];
vect_prior = {@IndicatifPrior, @IndicatifPrior, @IndicatifPrior};
[sample, gamma_sample, alpha_sample] = SamplesLogistic(@dfunLogistic, @MHMultVar_1, vect_prior, pos_guess, data_Train, vect_init, Time_step, gamma, T_0_sample,  T_sample, vect_alpha, vect_beta, opt_1);
acceptance_ratio = sum(alpha_sample(1,:))/length(alpha_sample);


for i = 1:nb_param
    mean_param(i) = mean(sample(i,:));
end
    
%Estimation with MH algorithm
vect_est = vect_init;
vect_est(pos_guess) = mean_param;
vect_temp = [vect_est(1) + mean(data_Train(1,:))/vect_est(3), vect_est(2), vect_est(3), vect_est(4)];
sol = ode45(@(t,y) dfunLogistic(t, y, vect_temp), tspan, mean(data_Train(1,:)), opt_1);
x_est = deval(sol, Time_step);

figure(num_fig)
std_Test = std(data_Train,0,2);
mean_Test = mean(data_Train, 2);
errorbar(Time_step, mean_Test, std_Test,'--o','MarkerSize',1,'LineWidth', 1,'MarkerEdgeColor','black','MarkerFaceColor','black');
hold on
plot(Time_step, x_est, '--', 'LineWidth', 2)
ylim([0 1.1*max(mean_Test)])
xlim([0 max(Time_step)])
legend('Simulation')
xlabel('Time category')
ylabel('Abundance (grams)')
title('Evolution abundance alone')
num_fig = num_fig + 1;

%%
%Determination of the distributions
%In this section, the goal is to find the underlying distribution
%generating the parameters mu_max and X_m. Then, these distributions will
%be used as prior for fittiong the parameters of our model.

mu_max_LN_p = lognfit(sample(1,:));
yield_LN_p = lognfit(sample(2,:));


%Histogram with log-normal fit for mu_{max}
figure(num_fig)
histogram(sample(1,:), 'Normalization', 'pdf');
hold on
fplot(@(x) lognpdf(x, mu_max_LN_p(1), mu_max_LN_p(2)), [0.95*min(sample(1,:)) 1.05*max(sample(1,:))], 'LineWidth', 2)%Log-Normal density
title('Log-Normal fit for mu_{max}')
num_fig = num_fig + 1;

%Histogram with log-normal fit yield
figure(num_fig)
histogram(sample(2,:), 'Normalization', 'pdf');
hold on
fplot(@(x) lognpdf(x, yield_LN_p(1), yield_LN_p(2)), [0.95*min(sample(2,:)) 1.05*max(sample(2,:))], 'LineWidth', 2)%Log-Normal density
title('Log-Normal fit for yield')
num_fig = num_fig + 1;

mean_y_0 = mean(data_Train(1,:));
std_y_0 = std(data_Train(1,:));

%% Second estimation

T_sample = 20000;%50000; %Length of the samples
T_0_sample = 200; %Remove first values, "calibration"
pos_guess = [2 3 4]; %[2 3 4];
nb_param = length(pos_guess);
k2_guess = exp(mu_max_LN_p(1) + mu_max_LN_p(2)^2/2) + normrnd(0, 0);%Mean of the log-normal distribution
k3_guess = exp(yield_LN_p(1) + yield_LN_p(2)^2/2) + normrnd(0, 0);%Mean of the log-normal distribution
k1_guess = (k2_guess + k3_guess)/mean_R_0;


mean_param = zeros(1,nb_param);
vect_init = [R_0_temp; k1_guess; k2_guess; k3_guess];
% Log-normal for yield and mu_max
vect_alpha = [1; 0; mu_max_LN_p(1) + normrnd(0, 0.1); yield_LN_p(1) + normrnd(0, 0.1)];
vect_beta = [1; 10^(6); 150*mu_max_LN_p(2); 150*yield_LN_p(2)];%[1; 10^(6); 50*mu_max_LN_p(2); 50*yield_LN_p(2)];
vect_prior = {@IndicatifPrior, @IndicatifPrior, @LNPrior, @LNPriork3};
sample = [k1_guess*ones(1, T_sample - T_0_sample); k2_guess*ones(1, T_sample - T_0_sample); k3_guess*ones(1, T_sample - T_0_sample)];
[sample((pos_guess - 1), :), gamma_sample_2, alpha_sample_2] = SamplesLogistic(@dfunTot, @MHMultVar_2, vect_prior, pos_guess, data_Train, vect_init, Time_step, gamma, T_0_sample,  T_sample, vect_alpha, vect_beta, opt_2);
acceptance_ratio_2 = sum(alpha_sample_2(1,:))/length(alpha_sample_2);

for i = 1:length(pos_guess)
    mean_param(i) = mean(sample(i,:));
end

%Estimation with MH algorithm
sol = ode45(@(t,x) dfunTot(t, x, mean(sample(1,:)), mean(sample(2,:)), mean(sample(3,:))), tspan, [mean(data_Train(1,:)) 0 0 mean_R_0], opt_2);
x_est = deval(sol, Time_step);
x_est = x_est(1,:);

figure(num_fig)
std_Test = std(data_Test,0,2);
mean_Test = mean(data_Test, 2);
errorbar(Time_step, mean_Test, std_Test,'--o','MarkerSize',1,'LineWidth', 1,'MarkerEdgeColor','black','MarkerFaceColor','black');
hold on
plot(Time_step, x_est, '--', 'LineWidth', 2)
set(gca, 'XTick', 0:10:100);
legend('Simulation')
xlabel('Time category')
ylabel('Abundance (grams)')
title('Evolution abundance alone')
num_fig = num_fig + 1;

%%
%Determination of the distributions
%In this section, the goal is to find the underlying distribution
%generating the parameters mu_max and X_m. Then, these distributions will
%be used to generate parameters for our simulations.

log_val_param = zeros(2,length(pos_guess));
for i = 1:length(pos_guess)
    log_val_param(:,i) = lognfit(sample(i,:));

    figure(num_fig)
    histogram(sample(i,:), 'Normalization', 'pdf');
    hold on
    fplot(@(x) lognpdf(x, log_val_param(1,i), log_val_param(2,i)), [0.95*min(sample(i,:)) 1.05*max(sample(i,:))], 'LineWidth', 2)%Log-Normal density
    title(strcat('Log-Normal fit for k', string(i)))
    num_fig = num_fig + 1;
end

CI_k1 = Conf_Interval(sample(1,randperm(T_sample - T_0_sample, 100)), 0.05);
CI_k2 = Conf_Interval(sample(2,randperm(T_sample - T_0_sample, 100)), 0.05);
CI_k3 = Conf_Interval(sample(3,randperm(T_sample - T_0_sample, 100)), 0.05); %Not commented if kappa_1 is estimated

if save_data == 1
    FolderName = strcat(cd, '/figures/');
    FigList = findobj(allchild(0), 'flat', 'Type', 'figure');
    for iFig = 1:length(FigList)
        FigHandle = FigList(iFig);
        FigName   = num2str(get(FigHandle, 'Number'));
        FigName = strcat(FolderName, FigName, CopyNum);
        set(0, 'CurrentFigure', FigHandle);
        saveas(FigHandle, FigName, 'pdf');
    end   
    LN_k1 = [log_val_param(1,1) log_val_param(2,1)];
    LN_k2 = [log_val_param(1,2) log_val_param(2,2)];
    LN_k3 = [log_val_param(1,3) log_val_param(2,3)];
    Abund_y_0 = [mean_y_0 std_y_0];
    save(strcat('./Data/',CopyNum, 'LN', 'Tablev2.mat'),'LN_k1', 'LN_k2', 'LN_k3'); 
    save(DataFile);
end

%% Functions 

function dz = dfunLogistic(t, x, vect)
x(x<=0) = 0;
dz = x*vect(2)*(1 - (x/(vect(3)*vect(1))));
end

function dz = dfunTot(t, x, kappa_1, kappa_2, kappa_3)
%kappa_3 = kappa_2/0.375 - kappa_2;
x(x<=0) = 0;
dx_tild = -kappa_1*x(1)*x(4) + (2*kappa_2 + kappa_3)*x(2);
dy = kappa_1*x(1)*x(4) - (kappa_2 + kappa_3)*x(2);
dc = kappa_3*x(2);
df = -kappa_1*x(1)*x(4);
dz = [dx_tild dy dc df]';
end

function [sample, gamma_sample, alpha_sample] = SamplesLogistic(dfun, funMH, vect_prior, pos_guess, x_obs, vect_init, t,  gamma, T_0_sample,  T_sample, vect_alpha, vect_beta, opt)
    %This function will deterime the variance of our distribution, either
    %we can choose it if we want a final distribution with more variability
    %or we can compute it based on the data (variance of the replicates)
    %and using and inverse gamma distribution. However, doing this could
    %lead to some errors if the values of the repliactes are very close.
    %This function uses the function MHMultVar -> Likelihood

    nb_guess = length(pos_guess);
    sample = vect_init(pos_guess).*ones(nb_guess,T_sample);
    gamma_sample = gamma(1)*ones(1,T_sample);
    alpha_sample = zeros(nb_guess,T_sample);
    theta = 0.1*vect_init;
    
    for j = 2:T_sample
        vect_init(pos_guess) = sample(:,j-1);
        [sol_temp, gamma, alpha, theta] = funMH(pos_guess, vect_init, dfun, vect_prior, gamma, x_obs, t, vect_alpha, vect_beta, j, theta, opt);
        sample(:,j) = sol_temp(pos_guess);
        gamma_sample(j) = gamma(1);
        alpha_sample(:,j) = alpha;
    end
    
    sample = sample(:,(T_0_sample+1):T_sample);
end

function [vect_init, gamma, alpha, theta] = MHMultVar_1(pos_guess, vect_init, dfun, vect_prior, gamma, x_obs, t, vect_alpha, vect_beta, j, theta, opt)
%This function performs the MH algorithm. In order to do this, we compute
%the likelihood assuming the data are normally distributed. This function
%uses the Likelihood function.
vect_new = vect_init;
alpha = ones(1,length(pos_guess));
nu = 2*j^(-3/4);%min(1, 2*j^(-3/4));
for i = 1:length(pos_guess)
    fun_prior = vect_prior{pos_guess(i)};
    vect_tild = vect_init;
    vect_tild(pos_guess(i)) = vect_init(pos_guess(i)) + normrnd(0, theta(pos_guess(i)));%normrnd(0, 0.1*vect_init(pos_guess(i)));%
    vect_tild(vect_tild<=0) = 10^(-4);
%     vect_tild(2:3)
    L_y = logLikelihood1(dfun, t, vect_tild(1), vect_tild(2), vect_tild(3), vect_tild(4), gamma, x_obs, opt);
    L_y = log(fun_prior(vect_tild(pos_guess(i)),vect_alpha(pos_guess(i)),vect_beta(pos_guess(i)), vect_tild(2))) + L_y;
    L_x = logLikelihood1(dfun, t, vect_init(1), vect_init(2), vect_init(3), vect_init(4), gamma, x_obs, opt);
    L_x = log(fun_prior(vect_init(pos_guess(i)),vect_alpha(pos_guess(i)),vect_beta(pos_guess(i)), vect_tild(2))) + L_x;
    ratio = min(exp(L_y - L_x),1);%If L_y >= L_x, take the new value "vect_tild". Otherwise (L_y < L_x) take this value ("vect_tild") with a probability ratio.
    rand = unifrnd(0,1);
    if rand > ratio
        vect_tild(pos_guess(i)) = vect_init(pos_guess(i));
        alpha(i) = 0; %Reject candidate
    end
    vect_new(pos_guess(i)) = vect_tild(pos_guess(i)); %The different parameters are assumed to be independent
    logtheta = log(theta(pos_guess(i))) + nu*(ratio - 0.234);
    theta(pos_guess(i)) = exp(logtheta);
end
vect_init = vect_new;
end

function [vect_init, gamma, alpha, theta] = MHMultVar_2(pos_guess, vect_init, dfun, vect_prior, gamma, x_obs, t, vect_alpha, vect_beta, j, theta, opt)
%This function performs the MH algorithm. In order to do this, we compute
%the likelihood assuming the data are normally distributed. This function
%uses the Likelihood function.
vect_new = vect_init;
alpha = ones(1, length(pos_guess));
nu = 2*j^(-3/4);%min(1, 2*j^(-3/4));
for i = 1:length(pos_guess) %Not the same prior for kappa_3 which is indicatif
    fun_prior = vect_prior{pos_guess(i)};
    vect_tild = vect_init;
    vect_tild(pos_guess(i)) = vect_init(pos_guess(i)) + normrnd(0, theta(pos_guess(i)));
    vect_tild(vect_tild<=0) = 10^(-4);
%     vect_tild(2)
%     vect_tild(3:4)
    L_y = log(fun_prior(vect_tild(pos_guess(i)), vect_alpha(pos_guess(i)), vect_beta(pos_guess(i)), vect_init(3)));
    if ~isinf(L_y)
        L_y = L_y + logLikelihood2(dfun, x_obs, gamma, t, vect_tild(1), vect_tild(2), vect_tild(3), vect_tild(4), opt);
    end
    L_x = logLikelihood2(dfun, x_obs, gamma, t, vect_init(1), vect_init(2), vect_init(3), vect_init(4), opt);
	L_x = log(fun_prior(vect_init(pos_guess(i)), vect_alpha(pos_guess(i)), vect_beta(pos_guess(i)), vect_init(3)))...
        + L_x;%Assume the parameters follows a log-normal distribution
    ratio = min(exp(L_y - L_x),1);%If L_y > L_x, take the new value "vect_tild". Otherwise take this value with a probability ratio.
	rand = unifrnd(0,1);
    if rand > ratio
        vect_tild(pos_guess(i)) = vect_init(pos_guess(i));
        alpha(i) = 0; %Reject candidate
    end
	vect_new(pos_guess(i)) = vect_tild(pos_guess(i));
    logtheta = log(theta(pos_guess(i))) + nu*(ratio - 0.234);
    theta(pos_guess(i)) = exp(logtheta);
end
vect_init = vect_new;
end

function lh = logLikelihood1(dfun, t, R_0, mu_max, yield, K_s, gamma, x_obs, opt)
%This function computes the likelihood of normally distributed errors.
%ODE logistic
%x_est = fun(t,lambda_est);
%Used when we only know the ODE
vect_init = [R_0, mu_max, yield, K_s];
vect_init(1) = vect_init(1) + mean(x_obs(1,:))/vect_init(3);
sol = ode45(@(t,y) dfun(t, y, vect_init), [min(t) max(t)], mean(x_obs(1,:)), opt);
x_est = deval(sol, t);
x_est = x_est';
lh = -1/2*sum(sum((x_est - x_obs).^2./(gamma'.^2),2));
end   

function lh = logLikelihood2(dfun, x_obs, gamma, Time_step, R_0, kappa_1, kappa_2, kappa_3, opt)
%This function computes the likelihood of normally distributed errors.
%Used when we only know the ODE
tspan = [min(Time_step), max(Time_step)]; %Time interval in hours
sol = ode45(@(t,x) dfun(t, x, kappa_1, kappa_2, kappa_3), tspan, [mean(x_obs(1,:)) 0 0 R_0], opt);
x_est = deval(sol, Time_step);
x_est = x_est(1,:);%sum(x_est(1:2,:));%
x_est = x_est';
lh = -1/2*sum(sum((x_est - x_obs).^2./(gamma'.^2),2));                      
end

function L = IndicatifPrior(x,a,b,y)
if x >= a && x<= b
    L = 1;
else
    L = 0;
end
end


function L = LNPriork3(x,a,b,y)
x = y/(y + x);
L = lognpdf(x, a, b);
end

function L = LNPrior(x,a,b,y)
L = lognpdf(x, a, b);
end


function CI = Conf_Interval(sample, alpha)
n = length(sample);
nu = n - 1; %Degree of freedom
z = tinv([alpha/2 (1-alpha/2)] , nu);
SE = std(sample)/sqrt(n);
CI = mean(sample) + z*SE;
end