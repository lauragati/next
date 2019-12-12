function [Aa, Ab, As] = matrices_A_intrate_smoothing2(param,hx)
% should be a smoother version that's automated and works directly with
% Mathematica, materials12.m (bottom)

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
% hand:
cxb = [sig, 1-bet,-sig*bet];  % correct (12 Dec 2019)
cxs = sig*[1 0 0 0]*(eye(nx)-bet*hx)^(-1);  % correct (12 Dec 2019)

cpa = [(1-alph)*bet, kapp*alph*bet,0];  % correct (12 Dec 2019)
cps = [0 0 1 0]*(eye(nx)-alph*bet*hx)^(-1);  % correct (12 Dec 2019)

cis = [0 1 0 rho];  % correct (12 Dec 2019)

% The below is copied from Mathematica, materials12.m (bottom)
w = (1+kapp.*psipi.*sig+psix.*sig);

% The g-coefficient matrices are the elements of the A matrices and come
% from the model solution in Mathematica, materials12.m (bottom)
g_pia = (cpa+cpa.*psix.*sig).*w.^(-1);
g_pib = cxb.*kapp.*w.^(-1);
g_pis = (cps+cxs.*kapp+(-1).*cis.*kapp.*sig+cps.*psix.*sig).*w.^(-1);

g_xa = (-1).*cpa.*psipi.*sig.*w.^(-1);  % correct (12 Dec 2019)
g_xb = cxb.*w.^(-1); % correct (12 Dec 2019)
g_xs = (cxs+(-1).*cis.*sig+(-1).*cps.*psipi.*sig).*w.^(-1);  % correct (12 Dec 2019)

g_ia = cpa.*psipi.*w.^(-1);
g_ib = cxb.*(kapp.*psipi+psix).*w.^(-1);
g_is = (cis+cps.*psipi+cxs.*kapp.*psipi+cxs.*psix).*w.^(-1);

% Now assemble everything you input and got from Mathematica
Aa = vertcat(g_pia, g_xa, g_ia);
Ab = vertcat(g_pib, g_xb, g_ib);
As = vertcat(g_pis, g_xs, g_is);

