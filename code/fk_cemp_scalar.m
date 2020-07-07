function k = fk_cemp_scalar(param,hx,Aa,Ab,As,a,b,SIG,kt_1)
% scalar criterion version of fk_CEMP
% 7 July 2020 

thetbar = param.thetbar;
gbar = param.gbar;

phi = [a,b];
[F,G] = FG(param,hx,Aa,Ab,As,a,b);

thet = abs(SIG^-1 * (phi - [F,G]));
% Take the row pertaining to differences in the inflation forecast (see Notes 7 July 2020)
thet = max(thet(1,:));

% % An analogue of Lutkepohl's chi2-statistic:
% taubar = (phi - [F,G])' * SIG^-1 *(phi - [F,G]);

I = thet <= thetbar;
k = I.*(kt_1+1)+(1-I).*gbar^(-1);
