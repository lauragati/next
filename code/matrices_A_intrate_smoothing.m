function [Ap_RE, As_RE, Aa_LR, Ab_LR, As_LR] = matrices_A_intrate_smoothing(param,hx)

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
g_pia = (1-kapp*sig*psi_pi/w)*[(1-alph)*bet, kapp*alph*bet, 0];
g_xa  =  -sig*psi_pi/w*[(1-alph)*bet, kapp*alph*bet, 0];
g_pib = kapp/w*[sig-sig*bet*psi_pi, 1-bet-sig*bet*psi_x, 0];
g_xb = 1/w*[sig-sig*bet*psi_pi, 1-bet-sig*bet*psi_x, 0];
g_pis = (1-kapp*sig*psi_pi/w)*[0 0 1 0]*(eye(nx)-alph*bet*hx)^(-1) - kapp*sig/w*[-1 1 0 rho]*(eye(nx)-bet*hx)^(-1);
g_xs = -sig*psi_pi/w*[0 0 1 0]*(eye(nx)-alph*bet*hx)^(-1) - sig/w*[-1 1 0 rho]*(eye(nx)-bet*hx)^(-1);

Aa_LR = vertcat(g_pia, g_xa, psi_pi*g_pia + psi_x*g_xa);
Ab_LR = vertcat(g_pib, g_xb, psi_pi*g_pib + psi_x*g_xb);
As_LR = vertcat(g_pis, g_xs, psi_pi*g_pis + psi_x*g_xs + [0 1 0 rho]);


% RE and EE models
Ap_RE = [bet + kapp*sig*(1-bet*psi_pi)/(w),kapp/w, 0; ...
     sig*(1-bet*psi_pi)/(w), 1/w, 0;...
    (psi_pi*(bet + kapp*sig*(1-bet*psi_pi)/(w)) +psi_x*sig*(1-bet*psi_pi)/(w)), psi_x/w+kapp*psi_pi/w, 0];

As_RE = [kapp*sig/w, -kapp*sig/w, 1-kapp*sig*psi_pi/w, 0; ...
    sig/w, -sig/w , -sig*psi_pi/w, 0; ...
    psi_x*(sig/w) + psi_pi*(kapp*sig/w), psi_x*(-sig/w) + psi_pi*(-kapp*sig/w)+1, psi_x*(-sig*psi_pi/w) + psi_pi*(1-kapp*sig*psi_pi/w), rho];