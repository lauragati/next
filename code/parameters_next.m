function [param, set] = parameters_next
param.bet  = 0.99;%0.99
param.sig  = 1; %1 IES=1 (log utility, consistent with long-run growth)
param.alph = 0.5; %0.5 (prob that firm stuck with price), set to correspond to an expected duration of 2 quarters. Collard: estimates of alpha between 0.66-0.75 (adjust prices every 3-4 months)
param.eta = 1/4; %1/4, inverse of Frisch. won't matter here. The Frisch: 1 or 2 (micro) or 4 (macro) e.g. Basu Lec 9, slide 16 Mac.
param.om = 1.25;% 1.25 Subsumes all things that go into elasticity of marginal cost to output (including Frisch) Woodford. Interest and Prices, p. 172 (Table 3.1).
param.thet = 10; %10. price elasticity of demand, Woodford, taken from Chari Kehoe McGrattan 2000.
param.zeta = (param.om + 1/param.sig)/(1+param.om*param.thet); % parameter of strategic complementarity. If < 1, strat comp in price setting. If >1, strt subs.
param.kapp = param.zeta * (1-param.alph*param.bet)/param.alph; % Woodford. Interest and Prices, p. 187.
param.psi_x  = 0; %0
param.psi_pi = 1.5; %1.5
param.w = 1+param.sig*param.psi_x +param.kapp*param.sig*param.psi_pi;
param.gbar    = 0.145; % 0.145 param_correct CEMP
param.thetbar = 1;%1 or 0.029; % param_correct CEMP
param.rho_r = 0; %0
param.rho_i = 0.6; % 0.6 CEMP: 0.877, not a perfect mapping, but the MC process, standing in for demand. Too much: set 0.7 to get monpol shock to increase i on impact
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
set.psi_x  = 0; 
set.psi_pi = 1.5; 
set.w = 1+set.sig*set.psi_x +set.kapp*set.sig*set.psi_pi;
set.gbar    = 0.145; 
set.thetbar = 0.029;
set.rho_r = 0;
set.rho_i = 0;
set.rho_u = 0; 
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

