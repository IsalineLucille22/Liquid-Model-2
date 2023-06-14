function [Tot_Gram_Evol, Abund_sim_PVE, Abund_sim_PPU, kappa_mat, mean_PVE_fin_props, Final_Abund_Sim_Pve, Final_Abund_Obs_Pve, W1_used, W2_used] = Main_Fun(Model_fun, mean_R_0, std_R_0, Mat_y_0, Time_step, Name_file, Name_Sheet, num_ratio, DataFile, mu_max_Pve, mean_LN_k2_Pve, mean_LN_k1_Pve, mean_LN_k3_Pve, mean_LN_k1_Ppu, mean_LN_k2_Ppu, mean_LN_k3_Ppu, yield_Ppu, Threshold_values, save_data)

    close all %Close opened figures

    opts_1 = odeset('RelTol',1e-9,'AbsTol',1e-10,'NonNegative',1:11); %Options for odeSolver

    yield_Pve = exp(mu_max_Pve(1))/(exp(mu_max_Pve(1)) + exp(mean_LN_k3_Pve(1)));
  
    %Output parameters
    tspan = [min(Time_step) max(Time_step)]; %Time in hours
    Num_Rep = 6; %Number of replicates

    %Observed final abundances %(Pve:Ppu) 100:1, 10:1, 1:1, 1:10, 1:100
    Obs_fin_Abund = readtable('./Data/DataMaximev3.xlsx', 'Sheet', 'FinalAbund'); 
    Obs_fin_Abund = [table2array(Obs_fin_Abund(1, 2:6)); table2array(Obs_fin_Abund(2, 2:6))];
   
    %Initialization of final results
    Abund_sim_PVE = zeros(Num_Rep, 1);
    Abund_sim_PPU = zeros(Num_Rep, 1);
    Fin_W1_used = zeros(Num_Rep, 1);
    Fin_W2_used = zeros(Num_Rep, 1);
    num_fig = 1;
    n_iter = Num_Rep; %Number of iterations. Number of simulated "replicates" for each initial ratios.
    Final_Abund_Obs_Pve = zeros(1,length(num_ratio));
    Final_Abund_Sim_Pve = zeros(1,length(num_ratio));
    for k = 1:length(num_ratio)
        Mat_y_0_temp = Mat_y_0(num_ratio(k),:);
        data = readtable('./Data/DataMaximev3.xlsx', 'Sheet', Name_Sheet{k});
        data_Evol = table2array(data);
        Threshold = Threshold_values(:, num_ratio(k));
        kappa_2 = [lognrnd(mean_LN_k2_Pve(1), mean_LN_k2_Pve(2),  1, n_iter);...0.363*ones(1, n_iter);...
            lognrnd(mean_LN_k2_Ppu(1), mean_LN_k2_Ppu(2),  1, n_iter)];
        kappa_3 = [kappa_2(1,:)/yield_Pve - kappa_2(1,:);...[lognrnd(mean_LN_k3_Pve(1), mean_LN_k3_Pve(2),  1, n_iter);...
            kappa_2(2,:)/yield_Ppu - kappa_2(2,:)];
%             lognrnd(mean_LN_k3_Ppu(1), mean_LN_k3_Ppu(2),  1, n_iter)];
        kappa_1_temp = [lognrnd(mean_LN_k1_Pve(1), mean_LN_k1_Pve(2),  1, n_iter);...
            lognrnd(mean_LN_k1_Ppu(1), mean_LN_k1_Ppu(2),  1, n_iter)];
        mean_R_0_temp = normrnd(mean_R_0, std_R_0);
        Tab_output = {};
        Tab_WastUsed = {};
        [Wastes_Evol_1, Wastes_Evol_2, Res_1_Evol, Res_2_Evol, X_1, X_2] = deal(zeros(length(Time_step),n_iter));
        for i = 1:n_iter
            mat_y_0 = [Mat_y_0_temp(1) Mat_y_0_temp(2) 0 0 0 0 mean_R_0_temp(1) mean_R_0_temp(2) 0 0 0];
            kappa_mat = [kappa_1_temp(1,i) kappa_2(1,i) kappa_3(1,i) kappa_1_temp(1,i); kappa_1_temp(2,i) kappa_2(2,i) kappa_3(2,i) kappa_1_temp(2,i); 3.374169e+03 0.09430373 0.2853004 3.374169e+03];
            sol = ode45(@(t, y) Model_fun(t, y, kappa_mat, Threshold), tspan,  mat_y_0, opts_1);
            z = deval(sol, Time_step);
            X_1(:,i) = z(1, :) + z(3, :); 
            X_2(:,i) = z(2, :) + z(4, :); 
            Tab_output{i} = [Time_step; X_1(:,i)'; X_2(:,i)'];
            Tab_WastUsed{i} = [z(9,:); z(10,:)]; %[Use of waste W1 by S2; Use of waste W2 by S1]
            Wastes_Evol_1(:,i) = z(5,:); %W_1 evolution
            Wastes_Evol_2(:,i) = z(6,:); %W_2 evolution
            Res_1_Evol(:,i) = z(7,:); %R evolution
            Res_2_Evol(:,i) = z(8,:); %R evolution
        end
        %Waste mean biomass evolution relative to time. Bar indicate the
        %standard deviation.
        figure(num_fig)
        num_fig = num_fig + 1;
        mean_Wastes_Evol_1 = mean(Wastes_Evol_1, 2); mean_Wastes_Evol_2 = mean(Wastes_Evol_2, 2);
        mean_X_Evol_1 = mean(X_1, 2); mean_X_Evol_2 = mean(X_2, 2);
        mean_Res_1_Evol = mean(Res_1_Evol, 2); mean_Res_2_Evol = mean(Res_2_Evol, 2);
        std_Wastes_Evol_1 = std(Wastes_Evol_1, [], 2); std_Wastes_Evol_2 = std(Wastes_Evol_2, [], 2);
        std_X_Evol_1 = std(X_1, [], 2); std_X_Evol_2 = std(X_2, [], 2);
        std_Res_1_Evol = std(Res_1_Evol, [], 2); std_Res_2_Evol = std(Res_2_Evol, [], 2);
        plot(Time_step, mean_X_Evol_1, 'b-', 'LineWidth', 1) 
        errorbar(Time_step, mean_X_Evol_1, std_X_Evol_1, 'Color', [0.4660 0.6740 0.1880],'LineWidth', 1) 
        hold on
        plot(Time_step, mean_X_Evol_2, 'r-', 'LineWidth', 1) 
        errorbar(Time_step, mean_X_Evol_2, std_X_Evol_2, 'Color', [0 0.4470 0.7410],'LineWidth', 1) 
        hold on
        errorbar(Time_step, mean_Wastes_Evol_1, std_Wastes_Evol_1, 'Color', [0 1 0],'LineWidth', 1) %Used by S_2
        hold on 
        errorbar(Time_step, mean_Wastes_Evol_2, std_Wastes_Evol_2, 'Color', [0.3010 0.7450 0.9330], 'LineWidth', 1) %Used by S_1
        hold on
        errorbar(Time_step, mean_Res_1_Evol, std_Res_1_Evol, 'Color', 	[1 0 0],'LineWidth', 1) %Evolution D-Mannitol
        hold on
        errorbar(Time_step, mean_Res_2_Evol, std_Res_2_Evol, 'Color', [0.8500 0.3250 0.0980], 'LineWidth', 1) %Evolution Putrescine
        ylim([0 5e-04])
        num_fig = num_fig + 1;

        %Randomly select nb_sim values among the n_iter values simulated
        ind_sim = randperm(n_iter);
        ind_sim = ind_sim(1:Num_Rep);
    

        [Tot_Gram_Evol, Tot_Gram_Evol_Pve, Tot_Gram_Evol_Ppu] = deal(zeros(length(Time_step), Num_Rep));
        [Props_sim_PVE, PVE_fin_abs, PPU_fin_abs, W1_tot_used, W2_tot_used]  = deal(zeros(1, Num_Rep));
        for i = 1:Num_Rep
            Temp = Tab_output{ind_sim(i)};
            Temp = Temp';
            Temp_W = Tab_WastUsed{ind_sim(i)};
            Temp_W = Temp_W';
            n_temp = length(Temp(:,1));
            Pve_Gram = Temp(:,2)/4.95e-04; %Conversion to obtain OD from Maximev3
            Ppu_Gram = Temp(:,3)/2.81e-04; %Conversion to obtain OD from Maximev3
            PVE_fin_abs(i) = Temp(n_temp,2);
            PPU_fin_abs(i) = Temp(n_temp,3);
            W1_tot_used(i) = Temp_W(n_temp,1);
            W2_tot_used(i) = Temp_W(n_temp,2);
            Tot_Gram_Evol(:,i) = Pve_Gram + Ppu_Gram; %Total value in OD (without noise). For fluorescence data, multiply by 0 the species that doesn't broadcast the fluorescence we are interested in.
            Tot_Gram_Evol_Pve(:,i) = Pve_Gram; %Only when we want to have the two OD evolution
            Tot_Gram_Evol_Ppu(:,i) = Ppu_Gram; %Only when we want to have the two OD evolution
            Props_sim_PVE(i) = Temp(n_temp,2)/(Temp(n_temp,3)+Temp(n_temp,2));
        end
        %Observed final abundances
        Props_obs_PVE = [Obs_fin_Abund(1,num_ratio(k))./(Obs_fin_Abund(1,num_ratio(k)) + Obs_fin_Abund(2,num_ratio(k))) NaN*ones(1, Num_Rep - 1)];
        PVE_obs = [Obs_fin_Abund(1,num_ratio(k)) NaN*ones(1, Num_Rep - 1)];
        Final_Abund_Obs_Pve(num_ratio(k)) = Obs_fin_Abund(1,num_ratio(k));
        Final_Abund_Sim_Pve(num_ratio(k)) = mean(PVE_fin_abs);
        W1_used = mean(W1_tot_used); W2_used = mean(W2_tot_used);
    
        % Figures generation
        % Boxplots stationary proportions
        figure(num_fig)
        boxplot([Props_sim_PVE; Props_obs_PVE]', 'Labels',{'Sim','Obs'})
        hold on 
        scatter(2*ones(size(Props_obs_PVE)).*(1+(rand(size(Props_obs_PVE))-0.5)/10),Props_obs_PVE,'r','filled')
        xlabel('Category')
        ylabel('Proportions')
        ylim([0.45 0.65])
        title('Boxplots final proportions Pve');
    
        % Boxplots stationary absolute biomasses
        num_fig = num_fig + 1;
        figure(num_fig)
        boxplot([PVE_fin_abs; PVE_obs]', 'Labels',{'Sim','Obs'})
        hold on 
        scatter(2*ones(size(PVE_obs)).*(1+(rand(size(PVE_obs))-0.5)/10),PVE_obs,'r','filled')
        xlabel('Category')
        ylabel('Proportions')
        title('Boxplots final abundances Pve');
    %     ylim([1.4e-04 2.5e-04])
        num_fig = num_fig + 1;
    
        % Optical density (OD) evolution relative to time

        figure(num_fig)
        col_index = 4:9;
%         col_index_2 = 14:19;%4:9, For OD and Pve Fluo. 14:19, For Ppu Fluo
        std_obs = std(data_Evol(:, col_index),0,2);
        mean_obs = mean(data_Evol(:, col_index), 2);
        errorbar(data_Evol(:, 1), mean_obs, std_obs,'m--o','MarkerSize',1,'LineWidth', 1,'MarkerEdgeColor','black','MarkerFaceColor','black');
        hold on
%         std_obs = std(data_Evol(:, col_index_2),0,2);
%         mean_obs = mean(data_Evol(:, col_index_2), 2);
%         errorbar(data_Evol(:, 1), mean_obs, std_obs,'c--o','MarkerSize',1,'LineWidth', 1,'MarkerEdgeColor','black','MarkerFaceColor','black');
%         hold on
        std_sim = std(Tot_Gram_Evol,0,2);
        mean_sim = mean(Tot_Gram_Evol, 2);
        errorbar(Time_step, mean_sim, std_sim,'b--o','MarkerSize',1,'LineWidth', 1,'MarkerEdgeColor','black','MarkerFaceColor','black');
        hold on
        std_sim = std(Tot_Gram_Evol_Pve,0,2);
        mean_sim = mean(Tot_Gram_Evol_Pve, 2);
        errorbar(Time_step, mean_sim, std_sim,'r--o','MarkerSize',1,'LineWidth', 1,'MarkerEdgeColor','black','MarkerFaceColor','black');
        hold on
        std_sim = std(Tot_Gram_Evol_Ppu,0,2);
        mean_sim = mean(Tot_Gram_Evol_Ppu, 2);
        errorbar(Time_step, mean_sim, std_sim,'g--o','MarkerSize',1,'LineWidth', 1,'MarkerEdgeColor','black','MarkerFaceColor','black');
        ylim([0 1.2])

        Abund_sim_PVE(1:Num_Rep) = PVE_fin_abs;
        Abund_sim_PPU(1:Num_Rep) = PPU_fin_abs;
        Fin_W1_used(1:Num_Rep) = W1_tot_used;
        Fin_W2_used(1:Num_Rep) = W2_tot_used;
        mean_PVE_fin_props = mean(Props_sim_PVE);
    
        %Wastes used
        num_fig = num_fig + 1;
        figure(num_fig)
        bar(flip([mean(Fin_W2_used) mean(Fin_W1_used)]), 'stacked')
        ylim([0 2e-04])
        num_fig = num_fig + 1;
    
        if save_data == 1
            FolderName = '/Users/isalinelucille-guex/Documents/Liquid-models/PutrescineDMan/Figures/Without cross feeding/'; %Change the path according to the localization of the Figure folder
            FigList = findobj(allchild(0), 'flat', 'Type', 'figure');
            for iFig = 1:length(FigList)
                FigHandle = FigList(iFig);
                FigName = num2str(get(FigHandle, 'Number'));
                FigName = strcat(FolderName, FigName, Name_file);
                set(0, 'CurrentFigure', FigHandle);
                saveas(FigHandle, FigName, 'pdf');
            end
            save(DataFile)
        end
    end
    % Comparison between Pve absolute stationary abundance observed vs
    % simulated for the different initial ratios.
    figure(num_fig)
    plot(num_ratio, Final_Abund_Obs_Pve, 'b--o')
    hold on
    plot(num_ratio, Final_Abund_Sim_Pve, 'r--o')
    ylim([1.5e-04 2.5e-04])
end