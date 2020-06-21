% Smooth anchoring function using functional approximation for the gain
% specified in levels, not changes
% 21 June 2020
function [k,g_pi,g_pibar] = fk_smooth_approx_univariate(alph,x,fe,fegrid_fine, g_fe)

xx = fe;
k1 = ndim_simplex_eval(x,xx,alph);

k=1/k1;

if nargin > 4 % only calc approx derivatives if gradient specified
%%%%%% approximate derivatives of the anchoring function, 13 June 2020

% find which gridpoint comes closest to the current point xx
offgrid_fe = fegrid_fine-xx(2);
closest_fe = abs(offgrid_fe) == min(abs(offgrid_fe));

% and take the derivatives at that gridpoint
g_pi = g_fe(closest_k, closest_fe);
g_pibar = -g_pi;

end
