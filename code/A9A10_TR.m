function [z, resids] = A9A10_TR(param,hx,fa,fb,s,seq)
% equations A9, A10 in materials24, evaluating residuals depending on
% number of input
% 2 April 2020

sig  = param.sig;
alph = param.alph;
bet  = param.bet;
kapp = param.kapp;
psi_pi = param.psi_pi;
psi_x  = param.psi_x;
nx = size(hx,1);

stuff1 = [sig,1-bet,-sig*bet]; % fa(3) I hope is rational
% stuff1 = [sig,1-bet,0];

stuff2 = sig*[1,0,0]*(eye(nx)-bet*hx)^(-1);
stuff3 = [(1-alph)*bet,kapp*alph*bet,0];
stuff4 = [0,0,1]*(eye(nx)-alph*bet*hx)^(-1);


if numel(seq)==1
    i=seq;
    resA9  = 0;
    resA10 = 0;
    x = -sig*i + stuff1*fb + stuff2*s;
    pi = kapp*x + stuff3*fa + stuff4*s;
    resTR = i - psi_pi*pi - psi_x*x;
elseif numel(seq)==2
    x = seq(1);
    i = seq(2);
    resA10 = 0;
    resA9 = -x -sig*i + stuff1*fb + stuff2*s;
    pi = kapp*x + stuff3*fa + stuff4*s;
    resTR = i - psi_pi*pi - psi_x*x;
elseif numel(seq) ==3
    pi = seq(1);
    x  = seq(2);
    i  = seq(3);
    resA9 = -x -sig*i + stuff1*fb + stuff2*s;
    resA10 =-pi +kapp*x + stuff3*fa + stuff4*s;
    resTR = i - psi_pi*pi - psi_x*x;
end

z =[pi;x;i];
resids = [resA9; resA10; resTR];

