function k = fk(zbar, b, s, kt_1, param, setp, Aa, Ab, As)
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
n = 3; % fortunately, this is the dimension of everything
P = eye(n).*[rho_r, rho_i, rho_u]';
SIG = eye(n).*[sig_r, sig_i, sig_u]';

M1 = eye(n) - Aa/(1-alph*bet)-Ab/(1-bet);
M2 = b - Aa*b*(eye(n)-alph*bet*P)*(-1)*P -Ab*b*(eye(n)-bet*P)^(-1)*P -As*P;
thet = abs(M1*zbar + M2*s);
I = thet <= thetbar*(sig_r + sig_i + sig_u);
k = I.*(kt_1+1)+(1-I).*gbar^(-1);

