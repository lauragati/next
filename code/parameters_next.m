function [param, set] = parameters_next
param.bet  = 0.98;
param.sig  = 0.5; % IES, mean across studies is 0.5
param.alph = 0.5;
param.kapp = (1-param.alph)*(1-param.alph*param.bet)/param.alph;
param.psi_x  = 1.5; 
param.psi_pi = 1.5; 
param.w = 1+param.sig*param.psi_x +param.kapp*param.sig*param.psi_pi;
param.gbar    = 0.145; % param_correct CEMP
param.thetbar = 20;%0.029; % param_correct CEMP
param.rho_r = 0.9;
param.rho_i = 0.9; 
param.rho_u = 0.9; 
param.sig_r = 0.1;%0.1; %?
param.sig_i = 0.359; % = sig_e from CEMP, standing in for the demand shock
param.sig_u = 0.277; % = sig_mu from CEMP, the cost-push shock

set.bet  = 0.98;
set.sig  = 0.5; 
set.alph = 0.5;
set.kapp =(1-set.alph)*(1-set.alph*set.bet)/set.alph;
set.psi_x  = 1.5; 
set.psi_pi = 1.5; 
set.w = 1+set.sig*set.psi_x +set.kapp*set.sig*set.psi_pi;
set.gbar    = 0.145; 
set.thetbar = 0.029;
set.rho_r = 0.9;
set.rho_i = 0.9;
set.rho_u = 0.9; 
set.sig_r = 0.1; 
set.sig_i = 0.359; 
set.sig_u = 0.277; 
set.ne = 3; % the number of shocks
set.nxnl = 3; % the number of nonlinear states (k + zbar, at 3 each -- > 3x2)
set.nxl = 3; % the number of linear states (z)
set.ny = 3;  % number of observables (also z for now)
set.nxnl_2 = 2; % the second dimension of the nonlinear states (2)

