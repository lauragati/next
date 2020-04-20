function k = fk_CEMP(param,hx,a,b,SIG,kt_1)
% accomodates vector learning
% output: 
% fk = k
% 26 Jan 2020
thetbar = param.thetbar;
gbar = param.gbar;

phi = [a,b];
[Aa, Ab, As] = matrices_A_25_true_baseline(param,hx);
[F,G] = FG(param,hx,Aa,Ab,As,a,b);

thet = max(max(abs(SIG^-1 * (phi - [F,G]))));

I = thet <= thetbar;
k = I.*(kt_1+1)+(1-I).*gbar^(-1);
