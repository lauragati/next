function [fa, fb] = fafb_anal_constant_free(param, set, a, b, s, hx)
% with 'free' I mean that it's independent of n
% 19 Oct 2019

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

nx = size(hx,1);

fa = a/(1-alph*bet) + b*(eye(nx)-alph*bet*hx)^(-1)*s;
fb = a/(1-bet) + b*(eye(nx)-bet*hx)^(-1)*s;