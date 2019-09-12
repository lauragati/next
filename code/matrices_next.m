function [B1, B2] = matrices_next(param, set)

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
C = eye(ne); % for now

g_pia = (1-kapp*sig*psi_pi/w)*[(1-alph)*bet, kapp*alph*bet, 0];
g_xa  =  -sig*psi_pi/w*[(1-alph)*bet, kapp*alph*bet, 0];
g_pib = kapp/w*[sig-sig*bet*psi_pi, 1-bet-sig*bet*psi_x, 0];
g_xb = 1/w*[sig-sig*bet*psi_pi, 1-bet-sig*bet*psi_x, 0];
g_pis = (1-kapp*sig*psi_pi/w)*[0 0 1]*(eye(ne)-alph*bet*P)^(-1) - kapp*sig/w*[-1 1 0]*(eye(ne)-bet*P)^(-1);
g_xs = -sig*psi_pi/w*[0 0 1]*(eye(ne)-alph*bet*P)^(-1) - sig/w*[-1 1 0]*(eye(ne)-bet*P)^(-1);

A1 = vertcat(g_pia, g_xa, psi_pi*g_pia + psi_x*g_xa);
A2 = vertcat(g_pib, g_xb, psi_pi*g_pib + psi_x*g_xb);
A3 = vertcat(g_pis, g_xs, psi_pi*g_pis + psi_x*g_xs + [0 1 0]);


B1 = A1/(1-alph*bet) + A2/(1-bet);
B2 = A1*C*(eye(ne) -alph*bet*P)^(-1) + A2*C*(eye(ne) -bet*P)^(-1) + A3;