function fk = fk_fun(k_t_1,pibar,gam,Gam,sigeta,gbar,thetbar)
disp('deprecated, check where you''re using fk_fun')
% I think sigeta = sige, just so you know...
I= abs((1-gam)*(Gam-1)*pibar) <= thetbar*sigeta;

fk = I.*(k_t_1+1) + (1-I).*gbar^(-1);