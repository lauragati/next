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
%%

true_fun  = @(x) x.^3;
true_fun2 = @(x) abs(x);
a=-1;
b=1;

X = linspace(a,b); % like the true sample
f = true_fun(X);

fhat1 = interp(f,2); % this one just stretches it by a factor of 2
N = 5; % number of nodes
Xq = linspace(a,b,N); % query points, i.e. the nodes themselves
fhat2 = interp1(X,f,Xq); % linear interpolation
fhat3 = interp1(X,f,Xq, 'spline'); % cubic spline interpolation
fhat4 = interp1(X,f,Xq, 'pchip'); % shape-preserving cubic (hopefully spline) interpolation



%% Plot
figure
plot(X,f); hold on
plot(Xq, fhat2)
plot(Xq, fhat3, '--')
plot(Xq, fhat4, 'o')
legend('Approximand', 'Linear interpolation', 'Cubic spline', 'Shape-preserving cubic')