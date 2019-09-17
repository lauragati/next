function [fa, fb] = fafb_anal_constant(param, set, a, s)
% 17 Sept 2019

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
n = 3;
P = eye(n).*[rho_r, rho_i, rho_u]';
SIG = eye(n).*[sig_r, sig_i, sig_u]';


% without P
fa = a/(1-alph*bet) + (eye(n)-alph*bet*P)^(-1)*s;
fb = a/(1-bet) + (eye(n)-bet*P)^(-1)*s;