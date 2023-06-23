function dx = fun_Inhibition(t, z, kappa, threshold)
%z(1), dx_1: biomass bacterial species 1. z(2), dx_2: biomass bacterial
%species 2. z(3), dy_1: biomass complex. z(4), dy_2: biomass complex.
%z(5), dw_1 biomass waste produced by species 1 that can be used by species
%2. z(6), dw_2: biomass waste produced by species 2 that can be used by
%species 1. z(7), dr: biomass resource. z(10) df: product release by S_1,
%S_2 when killed by W_1 or W_2.
% Model: S_i + R -> P_i -> 2S_i, P_i <-> S_i + W_i, S_j + W_i -> F.
z(z<=0) = 0;
k_1 = 4.526359e+01; k_2 = 4.526359e+01; %Bilateral Hill factor k = 3.567477e+01;
T_W_2 = 1/(1 + (threshold(1)/z(6))^k_1);%(z(6) - threshold(1) > 0); %Threshold function use of W_2 by S_1
T_W_1 = 1/(1 + (threshold(2)/z(5))^k_2);%(z(5) - threshold(2) > 0); %Threshold function use of W_1 by S_2
dx_1 = (2*kappa(1,2) + kappa(1,3))*z(3) - kappa(1,1)*z(1)*z(7) - T_W_2*kappa(1,4)*z(1)*z(6);
dx_2 = (2*kappa(2,2) + kappa(2,3))*z(4) - kappa(2,1)*z(2)*z(7) - T_W_1*kappa(2,4)*z(2)*z(5);
dy_1 = -(kappa(1,2) + kappa(1,3))*z(3) + kappa(1,1)*z(1)*z(7);
dy_2 = -(kappa(2,2) + kappa(2,3))*z(4) + kappa(2,1)*z(2)*z(7);
dw_1 = kappa(1,3)*z(3) - T_W_1*kappa(2,4)*z(2)*z(5);
dw_2 = kappa(2,3)*z(4) - T_W_2*kappa(1,4)*z(1)*z(6);
dr = -kappa(1,1)*z(1)*z(7) - kappa(2,1)*z(2)*z(7);
du_1 = T_W_1*kappa(2,4)*z(2)*z(5); %Used by PPU
du_2 = T_W_2*kappa(1,4)*z(1)*z(6); %Used by PVE
df =  T_W_2*kappa(1,4)*z(1)*z(6) + T_W_1*kappa(2,4)*z(2)*z(5);
dx = [dx_1; dx_2; dy_1; dy_2; dw_1; dw_2; dr; du_1; du_2; df];

%% Description
% Model: S_i + R -> P_i -> 2S_i, P_i <-> S_i + W_i, S_j + W_i -> F. 
% Model with different wastes produced by each species. The waste produced by the other species inhibits the development of the species consuming it when its concentration is above a certain threshold. 
% The threshold depends on the species and the detection function takes indicative form:
% max(z(j) - threshold(i), 0)/(z(j) - threshold(i))