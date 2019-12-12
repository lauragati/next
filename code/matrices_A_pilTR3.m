function [Aa, Ab, As] = matrices_A_pilTR3(param,hx)
% There's no "2", reflecting the fact that this uses the structure of
% matrices_A_intrate_smoothing3.m
% Mathematica, materials12.m (very bottom, section after "Restart 12 Dec 2019")

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

% transform these two to be compatible with Mathematica
psipi = psi_pi;
psix = psi_x;

% The c-coefficient matrices are model-specific and need to be specified by
% hand: (These are from Notes, 12 Dec 2019)
cxb = [sig, 1-bet-sig*bet*psi_x, 0]; % new
cxs = -sig*[-1 1 0 rho psipi]*(eye(nx)-bet*hx)^(-1); % new

cpa = [(1-alph)*bet, kapp*alph*bet, 0]; % unchanged
cps = [0 0 1 0 0]*(eye(nx)-alph*bet*hx)^(-1); % new

cis = [0 1 0 rho psipi];  % new

% Given the c-matrices, the solution is calculated in Mathematica and
% copied from there:
Aa = [cpa;[0 0 0];[0 0 0]];
Ab = [cxb.*kapp.*(1+psix.*sig).^(-1);cxb.*(1+psix.*sig).^(-1);cxb.*psix.*(1+psix.*sig).^(-1)];
As = [cps+cxs.*kapp.*(1+psix.*sig).^(-1);cxs.*(1+psix.*sig).^(-1);(1+psix.*sig).^(-1).*(cis+cxs.*psix+cis.*psix.*sig)];