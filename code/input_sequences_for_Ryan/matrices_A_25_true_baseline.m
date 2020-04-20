% matrices_25_true_baseline.m
% solves for the A-matrices allowing to differentiate between
% whether you impose relation (*), that agents know the Taylor rule, for
% Ab. See also smat.m. The way I write the ALM here is
% z_t = Aa*fa_t + Ab*fb_t + As*s_t
% see Notes 10 April 2020.
% 11 April 2020
function [Aa_april, Ab_april, As_april] = matrices_A_25_true_baseline(param,hx)
kapp = param.kapp;
sig  = param.sig;
psi_pi = param.psi_pi;
psi_x  = param.psi_x;
w = param.w;

[s1, s2, s3, s4, s5] = smat(param,hx);

Aa_april = [(1+sig*psi_x)/w;-sig*psi_pi/w;psi_pi/w]*s3;
Ab_april = [kapp/w; 1/w; (psi_x+kapp*psi_pi)/w]*s1;
As_april = [(kapp/w*s2 +(1+sig*psi_x)/w *s4    - kapp*sig/w*s5) ; ...
    (1/w*s2 - sig*psi_pi/w*s4                  -sig/w*s5);...
    ((psi_x+kapp*psi_pi)/w*s2 + psi_pi/w*s4    +1/w*s5)];

