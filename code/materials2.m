% materials 2
% Start working on simulating the mixed CEMP - Preston model
% 4 Sep 2019
clearvars
close all

% Add all the relevant paths
current_dir = pwd;
cd ../ % go up 1 levels
basepath = pwd;
cd .. % go up another level to BC_Research
BC_researchpath = pwd;
toolpath = [BC_researchpath '/matlab_toolbox'];
export_figpath = [toolpath '/Export_Fig'];
figpath = [basepath '/figures'];
tablepath = [basepath '/tables'];
datapath = [basepath '/data'];

cd(current_dir)

addpath(basepath)
addpath(toolpath)
addpath(export_figpath)
addpath(figpath)
addpath(datapath)

% Variable stuff ---
print_figs    = 0;
do_old_plots  = 0;
if print_figs ==1
    output_table  = 1;
else
    output_table =0;
end

fs=20; % fontsize


%% Simulation
[param, setp] = parameters_next;
param = struct2array(param);
setp = struct2array(setp);
T = 500;
burnin = 0;

rng(0)
[B1, B2] = matrices_next(param, setp);
[xn_true, xl_true, Y, e] = simul_next(param,setp,T, burnin, B1, B2);
[nx_nl1, nx_nl2] = size(xn_true(:,:,1));
nx_l = size(xl_true,1);
ny = size(Y,1);


% Observables
figure
set(gcf,'color','w'); % sets white background color
set(gcf, 'Position', get(0, 'Screensize')); % sets the figure fullscreen

subplot(1,ny,1)
plot(Y(1,:),'linewidth', 2)
ax = gca; % current axes
ax.FontSize = fs;
grid on
grid minor
title('Inflation')
subplot(1,ny,2)
plot(Y(2,:),'linewidth', 2)
ax = gca; % current axes
ax.FontSize = fs;
grid on
grid minor
title('Output gap')
subplot(1,ny,3)
plot(Y(3,:),'linewidth', 2)
ax = gca; % current axes
ax.FontSize = fs;
grid on
grid minor
title('Nom. int. rate')

if print_figs ==1
    figname = ['materials2_observables']
    cd(figpath)
    export_fig(figname)
    cd(current_dir)
    close
end

% Nonlinear states
figure
set(gcf,'color','w'); % sets white background color
set(gcf, 'Position', get(0, 'Screensize')); % sets the figure fullscreen

nonlin_titles = {'Gain - \pi', 'LR-E - \pi', 'Gain - x', 'LR-E - x',  'Gain - i','LR-E - i'};

subplot_index = 0;
for i=1:nx_nl1
    for j=1:nx_nl2
        subplot_index = subplot_index +1;
        subplot(nx_nl1,nx_nl2,subplot_index)
        plot(squeeze(xn_true(i,j,:)),'linewidth', 2)
        ax = gca; % current axes
        ax.FontSize = fs;
        grid on
        grid minor
        title(nonlin_titles(subplot_index))       
    end
end

if print_figs ==1
    figname = ['materials2_nonlin_states']
    cd(figpath)
    export_fig(figname)
    cd(current_dir)
    close
end

