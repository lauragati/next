function [z] = A9A10(param,hx,fa,fb,s,seq)
% equations A9, A10 in materials24
% 6 April 2020

sig  = param.sig;
alph = param.alph;
bet  = param.bet;
kapp = param.kapp;
nx = size(hx,1);

stuff1 = [sig,1-bet,-sig*bet]; % fa(3) I hope is rational
% stuff1 = [sig,1-bet,0];

stuff2 = sig*[1,0,0]*(eye(nx)-bet*hx)^(-1);
stuff3 = [(1-alph)*bet,kapp*alph*bet,0];
stuff4 = [0,0,1]*(eye(nx)-alph*bet*hx)^(-1);


if numel(seq)==1
    i=seq;
    x = -sig*i + stuff1*fb + stuff2*s;
    pi = kapp*x + stuff3*fa + stuff4*s;
elseif numel(seq)==2
    x = seq(1);
    i = seq(2);
    pi = kapp*x + stuff3*fa + stuff4*s;
elseif numel(seq) ==3
    pi = seq(1);
    x  = seq(2);
    i  = seq(3);
end

z =[pi;x;i];


