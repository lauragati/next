function [fa, fb, fy] = fafbfy(at_1,bt_1,param,st)

n = size(st,1);
k = n+1;
alph = param.alph;
bet = param.bet;
C = param.C;
f = param.f;

fa = at_1*(1-alph*bet)^(-1) + bt_1*(eye(n)-alph*bet*C)^(-1)*st;
fb = at_1*(1-bet)^(-1) + bt_1*(eye(n)-bet*C)^(-1)*st;
fy = f'*(eye(n)-bet*C)^(-1)*st;
