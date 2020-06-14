% Smooth anchoring function using functional approximation
% For use in estimation of the anchoring function, or later
% 10 June 2020
function [k,g_pi,g_pibar] = fk_smooth_approx(alph,x,fe,kt_1,k1grid_fine, fegrid_fine, g_fe)

xx = [1/kt_1; fe];
k1 = ndim_simplex_eval(x,xx,alph);
if k1<0
    warning('k1<0')
end

k=1/k1;

if nargin > 4 % only calc approx derivatives if gradient specified
%%%%%% approximate derivatives of the anchoring function, 13 June 2020

% % Construct finer grids that come closer to including xx
% k1grid = x{1};
% fegrid = x{2};
% ng_fine = length(k1grid)*10;
% 
% k1max = max(k1grid);
% k1min = min(k1grid);
% femax = max(fegrid);
% femin = min(fegrid);
% 
% k1grid_fine = linspace(k1min,k1max,ng_fine);
% fegrid_fine = linspace(femin,femax,ng_fine);
% [xxgrid, yygrid] = meshgrid(k1grid_fine,fegrid_fine);
% 
% % evaluate the approximation all over the fine grid
% k1_fine = ndim_simplex_eval(x,[xxgrid(:)';yygrid(:)'],alph);
% % evaluate gradient on the finer grid
% [g_k,g_fe] = gradient(reshape(k1_fine,ng_fine,ng_fine));

% find which gridpoint comes closest to the current point xx
offgrid_k = k1grid_fine-xx(1);
closest_k = abs(offgrid_k) == min(abs(offgrid_k));
offgrid_fe = fegrid_fine-xx(2);
closest_fe = abs(offgrid_fe) == min(abs(offgrid_fe));

% and take the derivatives at that gridpoint
g_pi = g_fe(closest_k, closest_fe);
g_pibar = -g_pi;

end


