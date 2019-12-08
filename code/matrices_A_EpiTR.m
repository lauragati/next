function [Aa_LR, Ab_LR, As_LR] = matrices_A_EpiTR(param,hx)
% Notes 7 Dec 2019
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
rho = param.rho;
nx= size(hx,1); % need to figure this out later

% LR model
% only the next four are different compared to the pi-in-TR models
g_pia = [(1-alph)*bet, kapp*alph*bet, 0];
g_xa  = [0, 0, 0];
g_xb = 1/(1+sig*psi_x)*[sig-sig*psi_pi, 1-bet, 0];
g_pib = kapp*g_xb;

g_pis = (1-kapp*sig*psi_pi/w)*[0 0 1 0]*(eye(nx)-alph*bet*hx)^(-1) - kapp*sig/w*[-1 1 0 rho]*(eye(nx)-bet*hx)^(-1);
g_xs = -sig*psi_pi/w*[0 0 1 0]*(eye(nx)-alph*bet*hx)^(-1) - sig/w*[-1 1 0 rho]*(eye(nx)-bet*hx)^(-1);

Aa_LR = vertcat(g_pia, g_xa, psi_pi*g_pia + psi_x*g_xa);
Ab_LR = vertcat(g_pib, g_xb, psi_pi*g_pib + psi_x*g_xb);
As_LR = vertcat(g_pis, g_xs, psi_pi*g_pis + psi_x*g_xs + [0 1 0 rho]);

