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
% the first is unchanged from the interest rate smoothing model, the second
% is changed
coeffa = [(1-alph)*bet, kapp*alph*bet, 0];
coeffb = [sig-sig*psi_pi, 1-bet-sig*bet*psi_x, 0];

% the next four are correct (verified in Mathematica, materials12.nb)
g_pia = coeffa;
g_xa  = 0*coeffa;
g_pib = kapp/(1+sig*psi_x)*coeffb;
g_xb = 1/(1+sig*psi_x)*coeffb;

% the next three are unchanged from the interest rate smoothing model
coefsa = [0 0 1 0]*(eye(nx)-alph*bet*hx)^(-1);
coefsb = [-1 1 0 rho]*(eye(nx)-bet*hx)^(-1);
coefss = [0 1 0 rho];

% the next two are correct (verified in Mathematica, materials12.nb)
g_pis = coefsa - kapp*sig/(1+sig*psi_x) *coefsb;
g_xs = -sig/(1+sig*psi_x)* coefsb;

Aa_LR = vertcat(g_pia, g_xa, psi_pi*g_pia + psi_x*g_xa);
Ab_LR = vertcat(g_pib, g_xb, psi_pi*g_pib + psi_x*g_xb);
As_LR = vertcat(g_pis, g_xs, psi_pi*g_pis + psi_x*g_xs + coefss);

