% Rewrite the general model w/o the Taylor-rule: given a sequence of i,
% what are x and pi? 
% See materials17.tex for a model summary and Notes 25 March 2020
% 25 March 2020
function [pi_x] = pi_x_given_i(param,hx,fa,fb,s,i)
sig  = param.sig;
alph = param.alph;
bet  = param.bet;
kapp = param.kapp;
nx = size(s,1);
pi_x = [0, 1; -kapp, 1]^(-1) * [-sig*i + [sig,1-bet,-sig*bet]*fb + sig*[1,0,0]*(eye(nx)-bet*hx)^(-1)*s; ...
    [(1-alph)*bet,kapp*alph*bet,0]*fa + [0,0,1]*(eye(nx)-bet*hx)^(-1)*s];