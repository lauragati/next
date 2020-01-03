function [Aa, Ab, As] = matrices_A2(param,hx)
% Mathematica, materials12d.m 
% 3 jan 2020

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

% The d-coefficient matrices are model-specific and need to be specified by
% hand % CONT HERE
d1 = [(1-alph)*bet, kapp*alph*bet, 0];
d2 = [0 0 1 0]*(eye(nx)-alph*bet*hx)^(-1);
d6 = [0 1 0];  

% Given the c-matrices, the solution is calculated in Mathematica and
% copied from there:
Aa = [cpa.*(1+psix.*sig).*(1+kapp.*psipi.*sig+psix.*sig).^(-1);(-1).* ...
 cpa.*psipi.*sig.*(1+(kapp.*psipi+psix).*sig).^(-1);cpa.*psipi.*(1+ ...
  (kapp.*psipi+psix).*sig).^(-1)];
Ab = [cxb.*kapp.*(1+(kapp.*psipi+psix).*sig).^(-1);cxb.*(1+(kapp.* ...
psipi+psix).*sig).^(-1);cxb.*(kapp.*psipi+psix).*(1+(kapp.*psipi+ ...
  psix).*sig).^(-1)];
As = [(1+kapp.*psipi.*sig+psix.*sig).^(-1).*(cps+cxs.*kapp+cps.*psix.* ...
  sig);(cxs+(-1).*cps.*psipi.*sig).*(1+kapp.*psipi.*sig+psix.*sig) ...
  .^(-1);(1+(kapp.*psipi+psix).*sig).^(-1).*(cis+cps.*psipi+cxs.*( ...
  kapp.*psipi+psix)+cis.*(kapp.*psipi+psix).*sig)];

