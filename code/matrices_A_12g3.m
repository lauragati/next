function [Aa, Ab, As] = matrices_A_12g3(param,hx)
% Mathematica, materials12g3.nb
% 9 jan 2020
% A-matrices for il-model, with "optimal forecaster" info assumption.
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
Aa = [bet+(-1).*alph.*bet,alph.*bet.*kapp,0;0,0,0;(-1).*((-1)+alph).* ...
  bet.*psipi,alph.*bet.*kapp.*psipi,0];
Ab = [kapp.*sig,kapp+(-1).*bet.*kapp,0;sig,1+(-1).*bet,0;(kapp.*psipi+ ...
  psix).*sig,(-1).*((-1)+bet).*(kapp.*psipi+psix),0];
As = [ia(3,1)+kapp.*sig.*(ib(1,1)+(-1).*bet.^(-1).*ib(4,1)),ia(3,2)+ ...
  kapp.*sig.*(ib(1,2)+(-1).*bet.^(-1).*ib(4,2)),ia(3,3)+kapp.*sig.*( ...
  ib(1,3)+(-1).*bet.^(-1).*ib(4,3)),ia(3,4)+bet.^(-1).*kapp.*sig.*( ...
  1+bet.*ib(1,4)+(-1).*ib(4,4));sig.*(ib(1,1)+(-1).*bet.^(-1).*ib(4, ...
  1)),sig.*(ib(1,2)+(-1).*bet.^(-1).*ib(4,2)),sig.*(ib(1,3)+(-1).* ...
  bet.^(-1).*ib(4,3)),bet.^(-1).*sig.*(1+bet.*ib(1,4)+(-1).*ib(4,4)) ...
  ;psix.*sig.*ib(1,1)+psipi.*(ia(3,1)+kapp.*sig.*ib(1,1))+(-1).* ...
  bet.^(-1).*(kapp.*psipi+psix).*sig.*ib(4,1),1+psix.*sig.*ib(1,2)+ ...
  psipi.*(ia(3,2)+kapp.*sig.*ib(1,2))+(-1).*bet.^(-1).*(kapp.*psipi+ ...
  psix).*sig.*ib(4,2),psix.*sig.*ib(1,3)+psipi.*(ia(3,3)+kapp.*sig.* ...
  ib(1,3))+(-1).*bet.^(-1).*(kapp.*psipi+psix).*sig.*ib(4,3),rho+ ...
  psipi.*ia(3,4)+kapp.*psipi.*sig.*ib(1,4)+psix.*sig.*ib(1,4)+(-1).* ...
  bet.^(-1).*(kapp.*psipi+psix).*sig.*((-1)+ib(4,4))];

mathematica_code = "/Users/lauragati/Dropbox/BC_Research/next/code/materials12g3.nb";
this_code = mfilename;
mname = strrep(this_code,'matrices_A_','')
mmcaname = extractBetween(mathematica_code,"/Users/lauragati/Dropbox/BC_Research/next/code/materials",".nb")

if strcmp(mname,mmcaname) ~=1
    error('not using the correct A-matrices from Mathematica!')
end
  