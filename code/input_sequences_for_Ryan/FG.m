function [F,G] = FG(param,hx,Aa,Ab,As,a,b)
% See Notes 25 Jan 2020
alph = param.alph;
bet=param.bet;
nx=size(hx,1);
F = ( Aa*1/(1-alph*bet) + Ab*1/(1-bet) ) * a;
G = Aa*b*(eye(nx)-alph*bet*hx)^(-1) +Ab*b*(eye(nx)-bet*hx)^(-1) + As;