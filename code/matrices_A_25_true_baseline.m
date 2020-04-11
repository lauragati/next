% matrices_25_true_baseline.m
% solves for the A-matrices w/o Mathematica, suspecting that the A-matrices
% from materials13 may be flawed (in particular, As is the suspect)
% see Notes 10 April 2020.
% 11 April 2020
function [Aa_april, Ab_april, As_april] = matrices_A_25_true_baseline(param,hx)
kapp = param.kapp;
sig  = param.sig;
psi_pi = param.psi_pi;
psi_x  = param.psi_x;
w = param.w;

[stuff1, stuff2, stuff3, stuff4] = stuff(param,hx);


aleph = [0,1,sig; 1,-kapp,0; -psi_pi,-psi_x,1];
aleph1 = [kapp/w, (1+sig*psi_x)/w, -kapp*sig/w; 1/w, -sig*psi_pi/w, -sig/w; (psi_x+kapp*psi_pi)/w, psi_pi/w, 1/w];

inv(aleph) - aleph1

Aa_april = [(1+sig*psi_x)/w;-sig*psi_pi/w;psi_pi/w]*stuff3;
Ab_april = [kapp/w; 1/w; (psi_x+kapp*psi_pi)/w]*stuff1;
As_april = [(kapp/w*stuff2 +(1+sig*psi_x)/w *stuff4) ; ...
    (1/w*stuff2 - sig*psi_pi/w*stuff4);...
    ((psi_x+kapp*psi_pi)/w*stuff2 + psi_pi/w*stuff4)];

