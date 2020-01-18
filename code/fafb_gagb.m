function [fa, fb, ga, gb] = fafb_gagb(param, a, b, c, d, s, hx)
% 18 Jan 2020

bet = param.bet;  
alph = param.alph;

nx = size(hx,1);

ga = 1/(1-alph*bet)*c + (eye(nx)-alph*bet*d)^(-1)*d*s;
gb = 1/(1-bet)*c + (eye(nx)-bet*d)^(-1)*d*s;
fa = a/(1-alph*bet) + b*ga;
fb = a/(1-bet) + b*gb;
