% Smooth anchoring function for the case that only pi is learned
% 19 March 2020
function [k, g_pi] = fk_smooth_pi_only(param,fe, kt_1)
% need to figure out how to input c, d 
d=param.d;
c=param.c;
% k1 = kt_1^(-1) + c + d*(fe^2);
% k = k1^(-1);
% % k=k1;
% k = kt_1 + d*1/(fe^2);
% k = kt_1 + 1/(d*fe^2);
k = kt_1 + 1/((d*fe)^2);
g_pi = -2*(d*fe)^(-2)*fe^(-1);