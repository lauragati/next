function [fb_straight, fb_queer]=fbtest(H,a,b,param,s)
n = size(s,1);
k = n+1;
bet = param.bet;
C = param.C;
f = param.f;

fb_straight = a*(1-bet)^(-1) + b*(eye(n)-bet*C)^(-1)*s;

h_fcst = zeros(k,H);
for h=1:H
%     h_fcst(:,h) = bet^h * (a + b*C^h *s);
    h_fcst(:,h) = bet^(h-1) * (a + b*C^(h-1) *s);
end
fb_queer = sum(h_fcst,2);