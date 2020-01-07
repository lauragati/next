function [Aa, Ab, As] = matrices_A_12f3(param,hx)
% Mathematica, materials12f3.nb
% 7 jan 2020
% A-matrices for il-model.
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

% The only thing you need to specify by hand is 
ia = (eye(nx)-alph*bet*hx)^(-1);
ib = (eye(nx)-bet*hx)^(-1);

% Given the d-matrices, the solution is calculated in Mathematica and
% copied from there:
Aa = 
Ab = 
As = 