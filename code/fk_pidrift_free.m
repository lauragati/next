function k = fk_pidrift_free(zbar, b, s, kt_1, param, Aa, Ab, As,hx)
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
nx = size(hx,1);

M1 = eye(nx-1) - Aa/(1-alph*bet)-Ab/(1-bet);
M2 = b - Aa*b*(eye(nx)-alph*bet*hx)*(-1)*hx -Ab*b*(eye(nx)-bet*hx)^(-1)*hx -As*hx;
thet = abs(M1*zbar + M2*s);
thet = thet(1);
I = thet <= thetbar*(sig_r + sig_i + sig_u);
k = I.*(kt_1+1)+(1-I).*gbar^(-1);
