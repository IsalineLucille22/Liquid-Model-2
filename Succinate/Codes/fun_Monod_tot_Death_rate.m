function dx = fun_Monod_tot_Death_rate(t, z, kappa, threshold)
%z(1), dx_1: biomass bacterial species 1. z(2), dx_2: biomass bacterial
%species 2. z(3), dy_1: biomass complex. z(4), dy_2: biomass complex.
%z(5), dw_1 biomass waste produced by species 1 that can be used by species
%2. z(6), dw_2: biomass waste produced by species 2 that can be used by
%species 1. z(7), dr: biomass resource.
% Model: S_i + R -> P_i -> 2S_i, P_i <-> S_i + W_i, S_j + W_i -> P_j, S_i
% -> 0. The last reaction corresponds to cell death rate (alpha).
z(z<=0) = 0;
T_1 = (z(6) - threshold(1) > 0);
T_2 = (z(5) - threshold(2) > 0);
alpha = [-0.1, -0.1];
alpha(1) = 0;%((5 <= t) && (t <= 11))*alpha(1);
alpha(2) = ((5 <= t) && (t <= 11))*alpha(2);
t_1 = 0; t_2 = 0; %Lag time, if there are lag times then t_i > 0.
dx_1 = (2*kappa(1,2) + kappa(1,3))*z(3) - kappa(1,1)*z(1)*z(7) - T_1*kappa(1,4)*z(1)*z(6) + alpha(1)*z(1); 
dx_2 = (2*kappa(2,2) + kappa(2,3))*z(4) - kappa(2,1)*z(2)*z(7) - T_2*kappa(2,4)*z(2)*z(5) + alpha(2)*z(1);
dy_1 = -(kappa(1,2) + kappa(1,3))*z(3) + kappa(1,1)*z(1)*z(7) + T_1*kappa(1,4)*z(1)*z(6);
dy_2 = -(kappa(2,2) + kappa(2,3))*z(4) + kappa(2,1)*z(2)*z(7) + T_2*kappa(2,4)*z(2)*z(5);
dw_1 = kappa(1,3)*z(3) - T_2*kappa(2,4)*z(2)*z(5);
dw_2 = kappa(2,3)*z(4) - T_1*kappa(1,4)*z(1)*z(6);
dr = -(t >= t_1)*kappa(1,1)*z(1)*z(7) - (t >= t_2)*kappa(2,1)*z(2)*z(7);
du_1 = T_2*kappa(2,4)*z(2)*z(5); %Used by PPU
du_2 = T_1*kappa(1,4)*z(1)*z(6); %Used by PVE
dx = [(t >= t_1)*dx_1; (t >= t_2)*dx_2; (t >= t_1)*dy_1; (t >= t_2)*dy_2; (t >= t_1)*dw_1; (t >= t_2)*dw_2; dr; du_1; du_2; 0];

%% Description
% Model: S_i + R -> P_i -> 2S_i, P_i <-> S_i + W_i, S_j + W_i -> P_j. 
% Model with different wastes produced by each species. The waste produced by the other species can be used when its concentration is above a certain threshold. 
% The threshold depends on the species and the detection function takes indicative form:
% max(z(j) - threshold(i), 0)/(z(j) - threshold(i))