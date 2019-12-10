function [Aa_LR, Ab_LR, As_LR] = matrices_A_pilTR(param,hx)
% 10 Dec 2019
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

% of the next five, only coeffb is changed wrt the interest rate smoothing
% model, and the coefs* guys get an extra element to account for the new
% state variable (it's always zero tho b/c pil doesn't enter any sum)
coeffa = [(1-alph)*bet, kapp*alph*bet, 0];
coeffb = [sig-sig*bet^2*psi_pi, 1-bet-sig*bet*psi_x, 0];  % <-- changed
coefsa = [0 0 1 0 0]*(eye(nx)-alph*bet*hx)^(-1);
coefsb = [-1 1 0 rho 0]*(eye(nx)-bet*hx)^(-1);
coefss = [0 1 0 rho 0];
coefspil = [0 0 0 0 1];

new_w = 1/(1+sig*psi_x +kapp*sig*psi_pi*bet); % <-- changed

% LR model - these all change
g_pia = (1-sig*psi_pi*bet*kapp/new_w) *coeffa;
g_xa  = -(sig*psi_pi*bet/new_w) *coeffa;
g_pib = kapp/new_w *coeffb;
g_xb = 1/new_w *coeffb;
g_pis = (1-kapp*sig*psi_pi*bet/new_w) *coefsa - kapp*sig/new_w *coefsb -kapp*sig*psi_pi/new_w *coefspil; % <-this is where pil shows up 1
g_xs = -sig*psi_pi*bet/new_w *coefsa - sig/new_w *coefsb -sig*psi_pi/new_w *coefspil; % <-this is where pil shows up 2

% here the nom interest rate won't respond to today's inflation, so the
% next three terms have a different last element
Aa_LR = vertcat(g_pia, g_xa, psi_x*g_xa);
Ab_LR = vertcat(g_pib, g_xb, psi_x*g_xb);
As_LR = vertcat(g_pis, g_xs, psi_x*g_xs + coefss + psi_pi*coefspil);  % <-this is where pil shows up 3
