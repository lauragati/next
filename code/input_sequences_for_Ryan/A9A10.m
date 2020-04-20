function [z] = A9A10(param,hx,fa,fb,s,seq)
% equations A9, A10 (in materials24 onwards)
% 6 April 2020

sig  = param.sig;
kapp = param.kapp;

[s1, s2, s3, s4] = smat(param,hx);


if numel(seq)==1
    i=seq;
    x = -sig*i + s1*fb + s2*s;
    pi = kapp*x + s3*fa + s4*s;
elseif numel(seq)==2
    x = seq(1);
    i = seq(2);
    pi = kapp*x + s3*fa + s4*s;
elseif numel(seq) ==3
    pi = seq(1);
    x  = seq(2);
    i  = seq(3);
end

z =[pi;x;i];


