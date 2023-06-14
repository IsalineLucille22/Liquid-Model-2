function dx = fun_Unique_Waste(t, z, kappa, threshold)
%z(1), dx_1: biomass bacterial species 1. z(2), dx_2: biomass bacterial
%species 2. z(3), dy_1: biomass complex. z(4), dy_2: biomass complex.
%z(5), dw biomass waste produced by species z(6), dw_2: biomass waste produced by species 2, equals 0 in this case z(7), dr: biomass resource.
% Model: S_i + R -> P_i -> 2S_i, P_i <-> S_i + W, S_j + W -> P_j.
z(z<=0) = 0;
t_1 = 0; t_2 = 0; %Lag time, if there are lag times then t_i > 0.
k = 35;
threshold = 1.2e-04;
kappa(1,4) = 20*kappa(1,1); kappa(2,4) = 0.1*kappa(2,1);
Hill_Coeff = (z(5) - threshold > 0); %1/(1 + (threshold/z(5))^k);%
dx_1 = (2*kappa(1,2) + kappa(1,3))*z(3) - kappa(1,1)*z(1)*z(7) - Hill_Coeff*kappa(1,4)*z(1)*z(5);
dx_2 = (2*kappa(2,2) + kappa(2,3))*z(4) - kappa(2,1)*z(2)*z(7) - Hill_Coeff*kappa(2,4)*z(2)*z(5);
dy_1 = -(kappa(1,2) + kappa(1,3))*z(3) + kappa(1,1)*z(1)*z(7) + Hill_Coeff*kappa(1,4)*z(1)*z(5);
dy_2 = -(kappa(2,2) + kappa(2,3))*z(4) + kappa(2,1)*z(2)*z(7) + Hill_Coeff*kappa(2,4)*z(2)*z(5);
dw_1 = kappa(1,3)*z(3) - Hill_Coeff*kappa(2,4)*z(2)*z(5)...
    + kappa(2,3)*z(4) - Hill_Coeff*kappa(1,4)*z(1)*z(5);
dw_2 = 0;
dr = -(t >= t_1)*kappa(1,1)*z(1)*z(7) - (t >= t_2)*kappa(2,1)*z(2)*z(7);
du_1 = Hill_Coeff*kappa(2,4)*z(2)*z(5); %Used by PPU
du_2 = Hill_Coeff*kappa(1,4)*z(1)*z(5); %Used by PVE
dx = [(t >= t_1)*dx_1; (t >= t_2)*dx_2; (t >= t_1)*dy_1; (t >= t_2)*dy_2; (t >= t_1)*dw_1; (t >= t_2)*dw_2; dr; du_1; du_2; 0];

%% Description
% Model: S_i + R -> P_i -> 2S_i, P_i <-> S_i + W_i, S_j + W_i -> P_j. 
% Model with one unique waste produced by each species. The waste produced can be used by both species with a different rate, when its concentration is above a certain threshold. 
% The threshold does not depend on the species and the detection function takes indicative form:
% max(z(j) - threshold(i), 0)/(z(j) - threshold(i))