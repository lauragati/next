% matrices_A_smat.m
% solves for the A-matrices allowing you to specify flexibly whether the
% Taylor rule is known and whether there is a monpol shock.
% Intended for the CEMP criterion and verified to work in
% clean_up_sim_learn.m and command_IRFs_approx_pretty.m
% 7 July 2020
function [Aa, Ab, As] = matrices_A_smat(param,hx, knowTR, mpshock)
kapp = param.kapp;
sig  = param.sig;
psi_pi = param.psi_pi;
psi_x  = param.psi_x;
w = param.w;

[s1, s2, s3, s4, s5] = smat(param,hx,knowTR,mpshock);


% this stuff comes from 11 April (Notes 9, p 27)
Aa = [(1+sig*psi_x)/w;-sig*psi_pi/w;psi_pi/w]*s3;
Ab = [kapp/w; 1/w; (psi_x+kapp*psi_pi)/w]*s1;
As = [(kapp/w*s2 +(1+sig*psi_x)/w *s4    - kapp*sig/w * s5) ; ...
    (1/w*s2 - sig*psi_pi/w*s4                  -sig/w * s5);...
    ((psi_x+kapp*psi_pi)/w*s2 + psi_pi/w*s4    +1/w * s5)];

