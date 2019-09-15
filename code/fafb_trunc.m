function [fa, fb] = fafb_trunc(param, set, phi, s, H)
% H is the horizon to which we're truncating here
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

a = phi(:,1);
b = phi(:,2:end);

Ea = zeros(n,H); % alph-bet discounted expectations
Eb = zeros(n,H); % bet discounted expectations
for h=1:H
    Ezh = Ez_h_general(param, set, phi, s, h);
    Ea(:,h) = (alph*bet)^(h-1) * Ezh;
    Eb(:,h) = (bet)^(h-1) * Ezh;
end
fa = sum(Ea,2);
fb = sum(Eb,2);