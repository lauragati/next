% play_interpolation_ND
% A playground for N-dimensional interpolation

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
%% Check the Curve Fitting Toolbox on the one-dimensional example

true_fun  = @(x) x.^3;
true_fun2 = @(x) abs(x);
a=-1;
b=1;

% Evaluation points (true sample)
Xe = linspace(a,b,20);
f = true_fun(Xe);
% f = true_fun2(Xe);
T = length(f);

% prepare nodes = query points Xq, and associated function values Yq
N = 5; % number of nodes
n = T/(N-1);
Xq = [Xe(1:n:end), Xe(end)];  % query points, i.e. the nodes themselves
Yq = [f(1:n:end), f(end)]; % function values at query points

% The original "spline" command
pps = spline(Xq,Yq);
fhat1 = ppval(pps,Xe);

ppc=csapi(Xq,Yq);
fhat2 = fnval(ppc,Xe);

figure
set(gcf,'color','w'); % sets white background color
set(gcf, 'Position', get(0, 'Screensize')); % sets the figure fullscreen
plot(Xe,f, 'linewidth',lw); hold on
plot(Xe, fhat1, 'linewidth',lw)
plot(Xe, fhat2, '--', 'linewidth',lw)
legend('Approximand', 'Cubic spline with "spline"', 'Cubic spline with "csapi"',...
    'location', 'northwest')
title('1D approximation')
ax = gca; % current axes
ax.FontSize = fs;
grid on
grid minor

%% Csapi is supposed to handle ND cases. So let's try a 2D case.

true_fun3  = @(x1,x2) x1+x2;
a=0;
b=1;
c=2;
d=3;

% Evaluation points and true f at evaluation points
T=10;
Xe = linspace(a,b,T);
Ye = linspace(c,d,T);
% using the notation that XX is the matrix version of the vector grid
[XXe,YYe] = ndgrid(Xe,Ye);
f = true_fun3(XXe,YYe);

% Interpolation points
N = 3; % number of nodes
n = T/(N-1);
Xq = [Xe(1:n:end), Xe(end)];
Yq = [Ye(1:n:end), Ye(end)];

[XXq, YYq]=ndgrid(Xq,Yq);
% evaluate function at query points
fq = true_fun3(XXq,YYq);
% evaluate spline at query points
pp = csapi({Xq,Yq},fq);
fhat3 = fnval(pp,{Xe,Ye});

figure
subplot(1,2,1)
mesh(XXe,YYe,f)
title('Approximand')
subplot(1,2,2)
mesh(Xe,Ye,fhat3)
title('Spline with "csapi"')
sgtitle('Indeed, look how beautifully csapi managed')