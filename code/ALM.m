% ALM.m
% Solves for the observables as a function of *s-matrices* and
% expectations
% see Notes 18 April 2020.
% 18 April 2020
function z = ALM(param,hx,fa,fb,s,knowTR,mpshock)

sig  = param.sig;
kapp = param.kapp;
psi_pi = param.psi_pi;
psi_x  = param.psi_x;

[s1, s2, s3, s4, s5] = smat(param,hx,knowTR,mpshock);

A = [0,1,sig; 1,-kapp,0; -psi_pi,-psi_x,1];

B = [s1*fb + s2*s; ...
    s3*fa + s4*s;...
    s5*s];

z=A\B;
