% computes the matrices I named "stuff" in front of fb, st, fa and st in
% equations A9 and A10 (e.g. materials25, appendix)
function [stuff1, stuff2, stuff3, stuff4, stuff5] = stuff(param,hx)
kapp = param.kapp;
bet  = param.bet;
alph = param.alph;
sig  = param.sig;
nx = length(hx);

stuff1 = [sig,1-bet,-sig*bet];
stuff2 = sig*[1,0,0]*(eye(nx)-bet*hx)^(-1);
stuff3 = [(1-alph)*bet,kapp*alph*bet,0];
stuff4 = [0,0,1]*(eye(nx)-alph*bet*hx)^(-1);
stuff5 = [0,1,0]; % if you wanna include a monpol shock.