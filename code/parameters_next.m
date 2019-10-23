function [param, set] = parameters_next
param.bet  = 0.98;%0.98
param.sig  = 0.5; % IES, mean across studies is 0.5
param.alph = 0.5; %0.5 (prob that firm stuck with price)
param.kapp = (1-param.alph)*(1-param.alph*param.bet)/param.alph;
param.psi_x  = 0; %1.5
param.psi_pi = 2; %1.5
param.w = 1+param.sig*param.psi_x +param.kapp*param.sig*param.psi_pi;
param.gbar    = 0.145; % 0.145 param_correct CEMP
param.thetbar = 1;%0.029; % param_correct CEMP
param.rho_r = 0; %0
param.rho_i = 0.877; %0.877, not a perfect mapping, but the MC process, standing in for demand
param.rho_u = 0; %0
param.rho = 0; % persistence of lag of interest rate
param.sig_r = 0.1;%0.1; %?
param.sig_i = 0.359; % 0.359 = sig_e from CEMP, standing in for the demand shock
param.sig_u = 0.277; % 0.277 = sig_mu from CEMP, the cost-push shock
param.kap =  0.80; % 0.8 allows the CUSUM test-statistic to be revised at a different rate than the estimate of the mean.  0 < kap < 1.
param.thettilde = 1.60; % 1.6 the new thetbar for CUSUM-test. I just set both to match CEMP's criterion approx.

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
set.rho = 0;
set.sig_r = 0.1; 
set.sig_i = 0.359; 
set.sig_u = 0.277; 
set.kap =  0.80; 
set.thettilde = 1.60; 
set.ne = 3; % the number of shocks
set.nxnl = 3; % the number of nonlinear states (k + zbar, at 3 each -- > 3x2)
set.nxl = 3; % the number of linear states (z)
set.ny = 3;  % number of observables (also z for now)
set.nxnl_2 = 2; % the second dimension of the nonlinear states (2)

