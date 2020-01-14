function [Aa, Ab, As] = matrices_A_12h2(param,hx)
% Mathematica, materials12h2.nb
% 10 jan 2020
% A-matrices for pil-model, with "myopic" info assumption.
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
Aa = [(1+(-1).*alph).*bet,alph.*bet.*kapp,0;0,0,0;0,0,0];
Ab = [(-1).*kapp.*sig.*((-1)+(-1).*psix.*sig).^(-1),(-1).*kapp.*((-1)+( ...
  -1).*psix.*sig).^(-1)+bet.*kapp.*((-1)+(-1).*psix.*sig).^(-1)+ ...
  bet.*kapp.*psix.*sig.*((-1)+(-1).*psix.*sig).^(-1),0;(-1).*sig.*(( ...
  -1)+(-1).*psix.*sig).^(-1),(-1).*((-1)+(-1).*psix.*sig).^(-1)+ ...
  bet.*((-1)+(-1).*psix.*sig).^(-1)+bet.*psix.*sig.*((-1)+(-1).* ...
  psix.*sig).^(-1),0;(-1).*psix.*sig.*((-1)+(-1).*psix.*sig).^(-1),( ...
  -1).*psix.*((-1)+(-1).*psix.*sig).^(-1)+bet.*psix.*((-1)+(-1).* ...
  psix.*sig).^(-1)+bet.*psix.^2.*sig.*((-1)+(-1).*psix.*sig).^(-1), ...
  0];
As = [ia(3,1)+(-1).*kapp.*sig.*((-1)+(-1).*psix.*sig).^(-1).*ib(1,1)+ ...
  kapp.*sig.*((-1)+(-1).*psix.*sig).^(-1).*ib(2,1)+kapp.*psipi.* ...
  sig.*((-1)+(-1).*psix.*sig).^(-1).*ib(4,1),ia(3,2)+(-1).*kapp.* ...
  sig.*((-1)+(-1).*psix.*sig).^(-1).*ib(1,2)+kapp.*sig.*((-1)+(-1).* ...
  psix.*sig).^(-1).*ib(2,2)+kapp.*psipi.*sig.*((-1)+(-1).*psix.*sig) ...
  .^(-1).*ib(4,2),ia(3,3)+(-1).*kapp.*sig.*((-1)+(-1).*psix.*sig).^( ...
  -1).*ib(1,3)+kapp.*sig.*((-1)+(-1).*psix.*sig).^(-1).*ib(2,3)+ ...
  kapp.*psipi.*sig.*((-1)+(-1).*psix.*sig).^(-1).*ib(4,3),ia(3,4)+( ...
  -1).*kapp.*sig.*((-1)+(-1).*psix.*sig).^(-1).*ib(1,4)+kapp.*sig.*( ...
  (-1)+(-1).*psix.*sig).^(-1).*ib(2,4)+kapp.*psipi.*sig.*((-1)+(-1) ...
  .*psix.*sig).^(-1).*ib(4,4);(-1).*sig.*((-1)+(-1).*psix.*sig).^( ...
  -1).*ib(1,1)+sig.*((-1)+(-1).*psix.*sig).^(-1).*ib(2,1)+psipi.* ...
  sig.*((-1)+(-1).*psix.*sig).^(-1).*ib(4,1),(-1).*sig.*((-1)+(-1).* ...
  psix.*sig).^(-1).*ib(1,2)+sig.*((-1)+(-1).*psix.*sig).^(-1).*ib(2, ...
  2)+psipi.*sig.*((-1)+(-1).*psix.*sig).^(-1).*ib(4,2),(-1).*sig.*(( ...
  -1)+(-1).*psix.*sig).^(-1).*ib(1,3)+sig.*((-1)+(-1).*psix.*sig).^( ...
  -1).*ib(2,3)+psipi.*sig.*((-1)+(-1).*psix.*sig).^(-1).*ib(4,3),( ...
  -1).*sig.*((-1)+(-1).*psix.*sig).^(-1).*ib(1,4)+sig.*((-1)+(-1).* ...
  psix.*sig).^(-1).*ib(2,4)+psipi.*sig.*((-1)+(-1).*psix.*sig).^(-1) ...
  .*ib(4,4);(-1).*psix.*sig.*((-1)+(-1).*psix.*sig).^(-1).*ib(1,1)+ ...
  psix.*sig.*((-1)+(-1).*psix.*sig).^(-1).*ib(2,1)+psipi.*psix.* ...
  sig.*((-1)+(-1).*psix.*sig).^(-1).*ib(4,1),(-1).*((-1)+(-1).* ...
  psix.*sig).^(-1)+(-1).*psix.*sig.*((-1)+(-1).*psix.*sig).^(-1)+( ...
  -1).*psix.*sig.*((-1)+(-1).*psix.*sig).^(-1).*ib(1,2)+psix.*sig.*( ...
  (-1)+(-1).*psix.*sig).^(-1).*ib(2,2)+psipi.*psix.*sig.*((-1)+(-1) ...
  .*psix.*sig).^(-1).*ib(4,2),(-1).*psix.*sig.*((-1)+(-1).*psix.* ...
  sig).^(-1).*ib(1,3)+psix.*sig.*((-1)+(-1).*psix.*sig).^(-1).*ib(2, ...
  3)+psipi.*psix.*sig.*((-1)+(-1).*psix.*sig).^(-1).*ib(4,3),(-1).* ...
  psipi.*((-1)+(-1).*psix.*sig).^(-1)+(-1).*psipi.*psix.*sig.*((-1)+ ...
  (-1).*psix.*sig).^(-1)+(-1).*psix.*sig.*((-1)+(-1).*psix.*sig).^( ...
  -1).*ib(1,4)+psix.*sig.*((-1)+(-1).*psix.*sig).^(-1).*ib(2,4)+ ...
  psipi.*psix.*sig.*((-1)+(-1).*psix.*sig).^(-1).*ib(4,4)];

mathematica_code = "/Users/lauragati/Dropbox/BC_Research/next/code/materials12h2.nb";
this_code = mfilename;
mname = strrep(this_code,'matrices_A_','')
mmcaname = extractBetween(mathematica_code,"/Users/lauragati/Dropbox/BC_Research/next/code/materials",".nb")

if strcmp(mname,mmcaname) ~=1
    error('not using the correct A-matrices from Mathematica!')
end
  