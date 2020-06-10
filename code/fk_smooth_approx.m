% Smooth anchoring function using functional approximation
% For use in estimation of the anchoring function, or later
% 10 June 2020
function k = fk_smooth_approx(alph,x,fe,kt_1)

xx = [1/kt_1; fe];
k1 = ndim_simplex_eval(x,xx,alph);




k=1/k1;
