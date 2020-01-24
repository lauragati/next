function [fa, fb] = fafb_materials14(param, a, b, s, hx)
% added a new hx-term in front of s
% 23 Jan 2020

bet = param.bet;  
alph = param.alph;

nx = size(hx,1);

fa = a/(1-alph*bet) + b*(eye(nx)-alph*bet*hx)^(-1)*hx*s;
fb = a/(1-bet) + b*(eye(nx)-bet*hx)^(-1)*hx*s;