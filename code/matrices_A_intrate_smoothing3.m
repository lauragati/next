function [Aa, Ab, As] = matrices_A_intrate_smoothing3(param,hx)
% a 2nd attempt at a smoother version that's automated and works directly with
% Mathematica, materials12.m (very bottom, section after "Restart 12 Dec 2019")
% Update 4 Jan 2020: I've checked that this is the correct MN method and
% indeed gives the same thing as the old A-matrices before I started using
% a general method. 

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
cxb = [sig-sig*bet*psi_pi, 1-bet-sig*bet*psi_x, 0];
cxs = -sig*[-1 1 0 rho]*(eye(nx)-bet*hx)^(-1);

cpa = [(1-alph)*bet, kapp*alph*bet, 0];
cps = [0 0 1 0]*(eye(nx)-alph*bet*hx)^(-1);

cis = [0 1 0 rho];  % correct (12 Dec 2019)

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

