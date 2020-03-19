% Smooth anchoring function for the case that only pi is learned
% 19 March 2020
function k = fk_smooth_pi_only(param,fe, kt_1)
% need to figure out how to input c, d - I'm cutting out c.
d=param.d;
k1 = kt_1^(-1) + d*fe^2;
k = k1^(-1);