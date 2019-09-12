function [fa, fb] = fafb_anal(param, set, phi, s)
% 12 Sept 2019

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
n = 3; % fortunately, this is the dimension of everything
P = eye(n).*[rho_r, rho_i, rho_u]';
SIG = eye(n).*[sig_r, sig_i, sig_u]';

fa = 1/(1-alph*bet)*phi + P*(eye(n)-alph*bet*P)^(-1)*s;
fb = 1/(1-bet)*phi + P*(eye(n)-bet*P)^(-1)*s;

% fa = 1/(1-alph*bet)*phi + (eye(n)-alph*bet*P)^(-1)*s;
% fb = 1/(1-bet)*phi + (eye(n)-bet*P)^(-1)*s;