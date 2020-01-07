function [Aa, Ab, As,  Ae] = matrices_A_12f1(param,hx)
% Mathematica, materials12f1.nb
% 7 jan 2020
% A-matrices for Epi-model.
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
Aa = [bet+(-1).*alph.*bet,alph.*bet.*kapp,0;0,0,0;0,0,0];
Ab = [(-1).*kapp.*((-1)+psipi).*sig.*(1+psix.*sig).^(-1),kapp.*((-1).* ...
  bet+(1+psix.*sig).^(-1)),0;(sig+(-1).*psipi.*sig).*(1+psix.*sig) ...
  .^(-1),(-1).*bet+(1+psix.*sig).^(-1),0;(-1).*((-1)+psipi).*psix.* ...
  sig.*(1+psix.*sig).^(-1),psix.*((-1).*bet+(1+psix.*sig).^(-1)),0]; ...
As = [ia(3,1)+kapp.*sig.*(1+psix.*sig).^(-1).*(ib(1,1)+(-1).*ib(2,1)), ...
  ia(3,2)+kapp.*sig.*(1+psix.*sig).^(-1).*(ib(1,2)+(-1).*ib(2,2)), ...
  ia(3,3)+kapp.*sig.*(1+psix.*sig).^(-1).*(ib(1,3)+(-1).*ib(2,3)); ...
  sig.*(1+psix.*sig).^(-1).*(ib(1,1)+(-1).*ib(2,1)),sig.*(1+psix.* ...
  sig).^(-1).*(ib(1,2)+(-1).*ib(2,2)),sig.*(1+psix.*sig).^(-1).*(ib( ...
  1,3)+(-1).*ib(2,3));psix.*sig.*(1+psix.*sig).^(-1).*(ib(1,1)+(-1) ...
  .*ib(2,1)),(1+psix.*sig).^(-1).*(1+psix.*sig.*(1+ib(1,2)+(-1).*ib( ...
  2,2))),psix.*sig.*(1+psix.*sig).^(-1).*(ib(1,3)+(-1).*ib(2,3))];
Ae = psipi;
