function [BM1, BM2] = bullard_mitra_matrices(param)
% [n,f,C,alph,bet,sig,kapp, psi_x,psi_pi,w_small] = parameters_preston;
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

BM1 = [bet + kapp*sig*(1-bet*psi_pi)/(w_small),kapp/w_small, 0; ...
    sig*(1-bet*psi_pi)/(w_small), 1/w_small, 0;...
    (psi_pi*(bet + kapp*sig*(1-bet*psi_pi)/(w_small)) +psi_x*sig*(1-bet*psi_pi)/(w_small)), psi_x/w_small+kapp*psi_pi/w_small, 0];

BM2 = [-kapp*sig/w_small, kapp/w_small; ...
    -sig/w_small, 1/w_small ; ...
    1-psi_x*sig/w_small-kapp*psi_pi*sig/w_small, psi_x/w_small+kapp*psi_pi/w_small];
