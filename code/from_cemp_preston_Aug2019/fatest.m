function [fa_straight, fa_queer]=fatest(H,a,b,param,s)
n = size(s,1);
k = n+1;
alph = param.alph;
bet = param.bet;
C = param.C;
f = param.f;

fa_straight = a*(1-alph*bet)^(-1) + b*(eye(n)-alph*bet*C)^(-1)*s;

h_fcst = zeros(k,H);
for h=1:H
%     h_fcst(:,h) = bet^h * (a + b*C^h *s);
    h_fcst(:,h) = (alph*bet)^(h-1) * (a + b*C^(h-1) *s);
end
fa_queer = sum(h_fcst,2);