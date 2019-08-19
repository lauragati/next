function [A1, A2, A3, A4] = preston_matrices(param)
n       = param.n;
f       = param.f;
C       = param.C; % IB is on top, RN is below
alph    = param.alph;
bet     = param.bet;
sig     = param.sig;
kapp    = param.kapp;
psi_x   = param.psi_x;
psi_pi  = param.psi_pi;
w_small = param.w_small;


A1 = [(1-kapp*sig*psi_pi/w_small)*[(1-alph)*bet, kapp*alph*bet, bet]; ...
    -sig*psi_pi/w_small*[(1-alph)*bet, kapp*alph*bet, bet];...
    (psi_pi*(1-kapp*sig*psi_pi/w_small) - psi_x*(sig*psi_pi/w_small))*[(1-alph)*bet, kapp*alph*bet, bet]];

A2 = [ kapp/w_small*[sig, 1-bet, -sig];...
    1/w_small*[sig, 1-bet, -sig];...
    (psi_pi*kapp/w_small + psi_x/w_small)*[sig, 1-bet, -sig]];

A3 = [ -kapp*sig/w_small*[1, -1];...
    -sig/w_small*[1, -1];...
    (-psi_pi*kapp*sig/w_small - psi_x*sig/w_small)*[1, -1]];
A4 = [0; -sig; 1];
