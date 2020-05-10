% play_interpolation
% A playground for interpolation, linear and spline, and shape-preserving

clearvars
close all
clc

% Add all the relevant paths and grab the codename
this_code = mfilename;
[current_dir, basepath, BC_researchpath,toolpath,export_figpath,figpath,tablepath,datapath] = add_paths;
todays_date = strrep(datestr(today), '-','_');

print_figs=0;
skip = 1;

% Plot configs
[fs, lw] = plot_configs;
%%

true_fun  = @(x) x.^3;
true_fun2 = @(x) abs(x);
a=-1;
b=1;

X = linspace(a,b,20); % like the true sample
% f = true_fun(X);
f = true_fun2(X);

T = length(f);

% "interp" command
fhat1 = interp(f,2); % this one just stretches it by a factor of 2


% prepare nodes = query points Xq, and associated function values Yq
N = 5; % number of nodes
n = T/(N-1);
Xq = [X(1:n:end), X(end)];  % query points, i.e. the nodes themselves
Yq = [f(1:n:end), f(end)]; % function values at query points

% the "interp1" command
fhat2 = interp1(Xq,Yq,X); % linear interpolation
fhat3 = interp1(Xq,Yq,X, 'spline'); % cubic spline interpolation
fhat4 = interp1(Xq,Yq,X, 'pchip'); % shape-preserving cubic hermite interpolating polynomial

% The "spline" command
% this gives you the coeffs, note that pp is a struct
pp = spline(Xq,Yq);
nodes = pp.breaks;
a = pp.coefs;
n_intervals =pp.pieces;
% this gives you the interpolant
fhat5 = spline(Xq,Yq,X);
% you can also do it using the coeffs explicity
fhat5_alt = ppval(pp,X);
fhat5-fhat5_alt; % bingo

% The "pchip" command - works exactly like the spline command
ppp= pchip(Xq,Yq);
fhat4_alt = ppval(ppp,X);
fhat4-fhat4_alt; % bingo

% Comparing pchip with SPLINE:
%     The function s(x) supplied by SPLINE is constructed in exactly the same way,
%     except that the slopes at the X(j) are chosen differently, namely to make 
%     even D^2s(x) continuous. This has the following effects.
%     SPLINE is smoother, i.e., D^2s(x) is continuous.
%     SPLINE is more accurate if the data are values of a smooth function.
%     pchip has no overshoots and less oscillation if the data are not smooth.
%     pchip is less expensive to set up.
%     The two are equally expensive to evaluate.


% let's compare with my thing
% set up objective spline with equation system
resids = obj_spline_natural(a,Xq,Yq);

% %Optimization Parameters
options = optimoptions('fsolve', 'TolFun', 1e-9, 'display', 'iter', 'MaxFunEvals', 10000);
tic
objh = @(coeffs) obj_spline_natural(coeffs,Xq,Yq);
[coeffs_opt,FVAL] = fsolve(objh,a, options);
toc

a-coeffs_opt % there's the problem
% although Matlab's spline uses the "not-a'knot" end conditions (specifies
% third derivates at endnodes) and not a natural spline, I have reason to
% suspect that my spline is not quite doing the right thing.

%% A multidimensional problem

true_fun3  = @(x1,x2) x1+x2;
a=0;
b=1;

X1 = linspace(a,b,20); 
X2 = linspace(a,b,20); 
% [X,Y] = ndgrid(X1,X2);
[X,Y] = meshgrid(X1,X2);

f = true_fun3(X,Y);


figure
mesh(X,Y,f)

return


%% Plot

figure
set(gcf,'color','w'); % sets white background color
set(gcf, 'Position', get(0, 'Screensize')); % sets the figure fullscreen
plot(X,f, 'linewidth',lw); hold on
plot(X, fhat2, 'linewidth',lw)
plot(X, fhat3, '--', 'linewidth',lw)
plot(X, fhat4, 'o', 'linewidth',lw)
plot(X, fhat5, 'x', 'linewidth',lw)
legend('Approximand', 'Linear interpolation with "interp1"', 'Cubic spline with "interp1"', 'Shape-preserving cubic with "interp1"', ...
    'Cubic spline with "spline"', ...
    'location', 'northwest')
ax = gca; % current axes
ax.FontSize = fs;
grid on
grid minor

%% FAZIT

% 1.) interp1('spline') = spline
% 2.) The logic in all these is to use a small number of nodes but pass a
% bigger number of data thru 'em
% 3.) Shape-preserving cubic spline thru the pchip option works well when
% there are kinks to avoid spline-like oscillations around 'em
% 4) there is a pchip command too> it's Piecewise Cubic Hermite Interpolating Polynomial