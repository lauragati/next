function k = fk_CEMP(param,hx,Aa,Ab,As,a,b,SIG,kt_1)
% reworked version of fk_pidrift_free
% accomodates vector learning
% output: 
% fk = k
% 26 Jan 2020
% sim_learnLH_clean was reworked 7 July 2020 to allow this code to have explicit assumptions on whether people know
% the Taylor rule and whether there is a monpol shock
thetbar = param.thetbar;
gbar = param.gbar;

phi = [a,b];
[F,G] = FG(param,hx,Aa,Ab,As,a,b);

thet = max(max(abs(SIG^-1 * (phi - [F,G]))));

% % An analogue of Lutkepohl's chi2-statistic:
% taubar = (phi - [F,G])' * SIG^-1 *(phi - [F,G]);

I = thet <= thetbar;
k = I.*(kt_1+1)+(1-I).*gbar^(-1);
