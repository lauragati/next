% Smooth anchoring function using functional approximation for the gain
% specified in levels, not changes. No inverse gains!
% 5 August 2020
function [k,g_pi,g_pibar] = fk_smooth_approx_univariate_uninvert(alph,x,fe,fegrid_fine, g_fe)

xx = fe;
k = ndim_simplex_eval(x,xx,alph);


if nargin > 4 % only calc approx derivatives if gradient specified
%%%%%% approximate derivatives of the anchoring function, 13 June 2020

% find which gridpoint comes closest to the current point xx
offgrid_fe = fegrid_fine-xx;
closest_fe = abs(offgrid_fe) == min(abs(offgrid_fe));

% and take the derivatives at that gridpoint
g_pi = g_fe(closest_fe);
g_pibar = -g_pi;

end
