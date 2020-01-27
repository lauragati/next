function [F,G] = FG_materials14(param,hx,Aa,Ab,As,a,b)
% FG with extra hx
% See Notes 25 Jan 2020
alph = param.alph;
bet=param.bet;
nx=size(hx,1);
F = ( Aa*1/(1-alph*bet) + Ab*1/(1-bet) ) * a;
G = Aa*b*(eye(nx)-alph*bet*hx)^(-1)*hx +Ab*b*(eye(nx)-bet*hx)^(-1)*hx + As;