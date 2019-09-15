function Ezh_t = Ez_h_general(param, set, phit_1, st, h)
% h-period ahead expectation (eq. 22 of materials3)
% 12 Sept 2019

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

a = phit_1(:,1);
b = phit_1(:,2:end);

% Ezh_t = phit_1*[1;P^(h-1)*st]; % this seems correct
% Ezh_t = phit_1*[1;P^(h)*st]; 
Ezh_t = a + b*P^(h-1)*st; % this seems to equal the top, correct one