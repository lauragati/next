function Epi_CB = Epi_CB(param, a, b, s, hx, Aa, Ab, As)
% 19 Jan 2020

bet = param.bet;  
alph = param.alph;

nx = size(hx,1);

Ez_CB = Aa * (a/(1-alph*bet) + b*(eye(nx)-alph*bet*hx)^(-1)*hx*s) ...
    + Ab * (a/(1-bet) + b*(eye(nx)-bet*hx)^(-1)*hx*s) ...
    + As * hx*s;

Epi_CB = Ez_CB(1,1);
