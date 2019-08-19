function [X1, X2] = lrexp2(fa,fb,fy,param)
alph = param.alph;
bet = param.bet;
kapp = param.kapp;
sig = param.sig;

X1 = [sig, 1-bet, -sig*bet]*fb + fy(2);
X2 = [(1-alph)*bet, kapp*alph*bet, 0]*fa;