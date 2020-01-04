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
d1 = sig*[1 -1 0]*(eye(nx)-bet*hx)^(-1);
d2 = [0 0 1 ]*(eye(nx)-alph*bet*hx)^(-1);
d6 = [0 1 0];  

% Given the d-matrices, the solution is calculated in Mathematica and
% copied from there:
Aa = [(-1).*((-1)+alph).*bet.*(1+psix.*sig).*(1+kapp.*psipi.*sig+psix.* ...
  sig).^(-1),alph.*bet.*kapp.*(1+psix.*sig).*(1+kapp.*psipi.*sig+ ...
  psix.*sig).^(-1),0;((-1)+alph).*bet.*psipi.*sig.*(1+kapp.*psipi.* ...
  sig+psix.*sig).^(-1),(-1).*alph.*bet.*kapp.*psipi.*sig.*(1+kapp.* ...
  psipi.*sig+psix.*sig).^(-1),0;(-1).*((-1)+alph).*bet.*psipi.*(1+ ...
  kapp.*psipi.*sig+psix.*sig).^(-1),alph.*bet.*kapp.*psipi.*(1+ ...
  kapp.*psipi.*sig+psix.*sig).^(-1),0];
Ab = [kapp.*(1+(-1).*bet.*psipi).*sig.*(1+kapp.*psipi.*sig+psix.*sig) ...
  .^(-1),(-1).*kapp.*(1+kapp.*psipi.*sig+psix.*sig).^(-1).*((-1)+ ...
  bet+bet.*psix.*sig),0;(sig+(-1).*bet.*psipi.*sig).*(1+kapp.* ...
  psipi.*sig+psix.*sig).^(-1),(-1).*(1+kapp.*psipi.*sig+psix.*sig) ...
  .^(-1).*((-1)+bet+bet.*psix.*sig),0;(-1).*((-1)+bet.*psipi).*( ...
  kapp.*psipi+psix).*sig.*(1+kapp.*psipi.*sig+psix.*sig).^(-1),(-1) ...
  .*(kapp.*psipi+psix).*(1+kapp.*psipi.*sig+psix.*sig).^(-1).*((-1)+ ...
  bet+bet.*psix.*sig),0];
As = [d1.*kapp.*(1+kapp.*psipi.*sig+psix.*sig).^(-1)+(-1).*d2.*((-1)+( ...
  -1).*psix.*sig).*(1+kapp.*psipi.*sig+psix.*sig).^(-1);d1.*(1+ ...
  kapp.*psipi.*sig+psix.*sig).^(-1)+(-1).*d2.*psipi.*sig.*(1+kapp.* ...
  psipi.*sig+psix.*sig).^(-1);d6+d2.*psipi.*(1+kapp.*psipi.*sig+ ...
  psix.*sig).^(-1)+d1.*(kapp.*psipi.*(1+kapp.*psipi.*sig+psix.*sig) ...
  .^(-1)+psix.*(1+kapp.*psipi.*sig+psix.*sig).^(-1))];
