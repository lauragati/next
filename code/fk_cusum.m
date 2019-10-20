function [k, om, thet] = fk_cusum(param,kt_1,omt_1, thett_1,f)
% alternative specification for theta_t as a CUSUM test (CEMP p.18 (p.19, Mac))
% 14 oct 2019
bet = param.bet;  
sig = param.sig;
alph = param.alph;
kapp = param.kapp;
psi_x = param.psi_x;
psi_pi = param.psi_pi;
w = param.w;
gbar = param.gbar;
thetbar = param.thetbar;
rho_r = param.rho_r;
rho_i = param.rho_i;
rho_u = param.rho_u;
sig_r = param.sig_r;
sig_i = param.sig_i;
sig_u = param.sig_u;
% n = 3; % fortunately, this is the dimension of everything
% P = eye(n).*[rho_r, rho_i, rho_u]';
% SIG = eye(n).*[sig_r, sig_i, sig_u]';

kap = 0.80; %0.035 allows the test-statistic to be revised at a different rate than the estimate of the mean. For now, I just set it to the mean value of CEMP's calibration.
thettilde = 1.60; % 0.75 the new thetbar, I set it to the mean of the calibration of CEMP.
om = omt_1 + kap*kt_1^(-1)*(f^2 - omt_1);
thet = thett_1 + kap*kt_1^(-1)*(f^2/om - thett_1);
I = thet <= thettilde;
k = I.*(kt_1+1)+(1-I).*gbar^(-1);

