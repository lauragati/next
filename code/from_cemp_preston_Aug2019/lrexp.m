function [X1, X2] = lrexp(H,a,b,param,s)
n = size(s,1);
k = n+1;
alph = param.alph;
bet = param.bet;
kapp = param.kapp;
sig = param.sig;
C = param.C;
f = param.f;

Ez = zeros(k,H);
Ebz = zeros(k,H);
Eabz = zeros(k,H);
Es = zeros(n,H);
for h=1:H
    % h-period ahead expectation of y
    Ez(:,h) = a + b*C^(h-1) *s;
    Ebz(:,h) = bet^(h-1) * Ez(:,h); % bet-discounted h-period ahead fcst
    Eabz(:,h) = (alph*bet)^(h-1) * Ez(:,h); % alph*bet-discounted h-period ahead fcst
    % h-period ahead expectation of rN, bet discounted
    Es(:,h) = bet^(h-1) * C^(h-1)*f'*s;
end

sumEb = sum(Ebz,2);
sumEab = sum(Eabz,2);
sumEs = sum(Es,2);

X1 = [sig, 1-bet, -sig*bet]*sumEb + sumEs(2);
X2 = [(1-alph)*bet, kapp*alph*bet, 0]*sumEab;