% aleph_gimel2.m
% aleph_gimel for a first-pass reaction function for i, r1.
% 13 April 2020
function z = aleph_gimel2(param,hx,fa,fb,s,k,pibar,xbar)

sig  = param.sig;
kapp = param.kapp;
psi_pi = param.psi_pi;
psi_x  = param.psi_x;
psi_k = param.psi_k;
psi_pibar  = param.psi_pibar;
psi_xbar  = param.psi_xbar;

[stuff1, stuff2, stuff3, stuff4,stuff5] = stuff(param,hx);

aleph = [0,1,sig; 1,-kapp,0; -psi_pi,-psi_x,1];

gimel = [stuff1*fb + stuff2*s; ...
    stuff3*fa + stuff4*s;...
    psi_k/k + psi_pibar*pibar + psi_xbar*xbar + stuff5*s];

z=aleph\gimel;