function [fa, fb] = fafb(param, set, zbar, s)

bet  = param(1);
sig  = param(2);
alph = param(3);
kapp = param(4);
psi_x  = param(5);
psi_pi = param(6);
w    = param(7);
gbar = param(8);
thetbar = param(9);
rho_r = param(10);
rho_i = param(11);
rho_u = param(12);
sig_r = param(13);
sig_i = param(14); 
sig_u = param(15); 
ne = set(16);
nxnl = set(17);
nxl = set(18);
ny = set(19);
P = eye(ne).*[rho_r, rho_i, rho_u]';

fa = 1/(1-alph*bet)*zbar + (eye(ne)-alph*bet*P)^(-1)*s;
fb = 1/(1-bet)*zbar + (eye(ne)-bet*P)^(-1)*s;