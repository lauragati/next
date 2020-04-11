% aleph_gimel.m
% Replaces matrices_25_true_baseline.m in not solving explicitly for the
% A-matrices, but directly for the observables as a function of *stuff* and
% expectations
% see Notes 11 April 2020.
% 11 April 2020
function z = aleph_gimel(param,hx,fa,fb,s)

sig  = param.sig;
kapp = param.kapp;
psi_pi = param.psi_pi;
psi_x  = param.psi_x;

[stuff1, stuff2, stuff3, stuff4,stuff5] = stuff(param,hx);

aleph = [0,1,sig; 1,-kapp,0; -psi_pi,-psi_x,1];

gimel = [stuff1*fb + stuff2*s; ...
    stuff3*fa + stuff4*s;...
    stuff5*s];

z=aleph\gimel;