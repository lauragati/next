function [z] = A9A10(param,hx,fa,fb,s,seq)
% equations A9, A10 in materials24
% 6 April 2020

sig  = param.sig;
kapp = param.kapp;

[stuff1, stuff2, stuff3, stuff4] = stuff(param,hx);


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


