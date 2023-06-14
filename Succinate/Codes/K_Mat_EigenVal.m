function [eig_val, eig_vect] = K_Mat_EigenVal(X_1, X_2, kappa_mat)
K_mat = [-(kappa_mat(2,4) + kappa_mat(2,3))*X_2 kappa_mat(2,4)*X_2;...
    kappa_mat(1,4)*X_1 -(kappa_mat(1,4) + kappa_mat(1,3))*X_1];
%[eig_vect, eig_val] = eig(K_mat);
K_star_mat = 1/2*(K_mat + K_mat');
[eig_vect, eig_val] = eig(K_star_mat);
end

%% Description
%Compute the eigen values and eigen vector of the matrix K for the
%simplified system without complex:
% S_i + R -> 2*S_i, S_i + R -> S_i + W_i, S_i + W_j -> 2*S_i, S_i + W_j ->
% S_i + W_i.