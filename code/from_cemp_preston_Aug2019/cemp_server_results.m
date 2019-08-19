% Load in output from server here to do smoother and generate figures and tables
% Add all the relevant paths
clearvars 

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
serverpath = [current_dir '/server'];

cd(current_dir)

addpath(basepath)
addpath(toolpath)
addpath(export_figpath)
addpath(figpath)
addpath(datapath)
addpath(serverpath)

% Variable stuff ---
print_figs    = 0;
do_old_plots  = 0;
if print_figs ==1
    output_table  = 1;
else
    output_table =0;
end
fs=20; % fontsize

load('data_cemp_timevarying.mat')
load('cemp_bayesian_estimation_outputs_server.mat')
load('cemp_bayesian_estimation_configs.mat') % N,short_chaint, long_chaint, accept_rate, duration

param_bay = mean(pchain);

if output_table==1
% Gather information for tables:
param_correct =  [2.472 0.029     0.145 0.128  0.891  0.877   0.084 0.359  0.277 0.042   0.021  0.073 0.049];
parameter_table = [param_correct; param_bay];
% pistar, thetbar, gbar, gamma, Gamma, rhophi, sige, sigmu, sigo1, sigo2, sigo3, sigo4,sigo5
columnLabels = {'$\pi^*$', '$\bar{\theta}$','$\bar{g}$','$\gamma$','$\Gamma$','$\rho$','$\sigma^2_e$','$\sigma^2_{\mu}$','$\sigma^2_{o_1}$','$\sigma^2_{o_2}$','$\sigma^2_{o_3}$','$\sigma^2_{o_4}$','$\sigma^2_{o_5}$'};
rowLabels = {'True', 'Posterior'};
matrix2latex_black(parameter_table, [tablepath '/materials8_paramtable_estim_cemp_server.tex'], 'rowLabels', rowLabels, 'columnLabels', columnLabels, ...
    'alignment', 'c', 'format', '%-6.3f', 'size', 'small');

estimation_configs = [N, short_chaint, long_chaint, accept_rate*100];
rowLabels = {[' Runtime: ', datestr(duration,'HH:MM:SS')]};
columnLabels = {'\# particles (N)', 'Short chain length', 'Long chain length', 'Acceptance rate (\%)'};
matrix2latex_black(estimation_configs, [tablepath '/materials8_estim_configs_cemp_server.tex'], 'rowLabels', rowLabels, 'columnLabels', columnLabels, ...
    'alignment', 'c', 'format', '%-6.0f', 'size', 'small');
end
%% Smoothed states
disp('--------- Marginalized particle smoother on server output ---------')
disp(['Started at hh:mm:ss: ', datestr(now,'HH:MM:SS')])

% Run the forward filter to get end-of-history swarms
[~, xn_swarm, xl_swarm, Ptt_seq, qt_seq, qtilde_seq] = MPF_test_missings(nanY,Y1,Y2,Y3,param_bay,N);

M = 100; % number of trajectories of smoothed states
tic
% dbstop in marginalized_smoother_test at 76 if t==470
[pibar_smooth, k_smooth, xi_smooth] = marginalized_smoother_test(nanY,param_bay,N,M,xn_swarm,xl_swarm,Ptt_seq,qtilde_seq);
toc

savefilename = 'smoothed_states_cemp_estim_server.mat';
save(savefilename, 'pibar_smooth','k_smooth', 'xi_smooth')

%% Plotting smoothed states
load('smoothed_states_cemp_estim_server.mat')
meanpibar = mean(pibar_smooth,2);
meank = mean(k_smooth,2);
meanxi = squeeze(mean(xi_smooth, 3));

figure
set(gcf,'color','w'); % sets white background color
set(gcf, 'Position', get(0, 'Screensize')); % sets the figure fullscreen
subplot(1,2,1)
plot(datenum_data,meanpibar, 'k', 'linewidth',2); hold on
ax = gca; % current axes
ax.FontSize = fs;
grid on
grid minor
legend('Long-term expectations, smoothed estimate')
datetick('x','yyyy-qq')

subplot(1,2,2)
plot(datenum_data, 1./meank, 'k', 'linewidth',2); hold on
ax = gca; % current axes
ax.FontSize = fs;
grid on
grid minor
legend('Inverse gain, smoothed estimate')
datetick('x','yyyy-qq')
figname = ['materials8_cemp_fig2_server'];
if print_figs ==1
    cd(figpath)
    export_fig(figname)
    cd(current_dir)
    close
end

figure
set(gcf,'color','w'); % sets white background color
set(gcf, 'Position', get(0, 'Screensize')); % sets the figure fullscreen
plot(datenum_data, meanxi(:,3), 'k', 'linewidth',2); hold on
plot(datenum_data, nanY(1,:), 'b', 'linewidth',2); hold on
ax = gca; % current axes
ax.FontSize = fs;
grid on
grid minor
legend('Inflation, smoothed estimate', 'CPI inflation')
datetick('x','yyyy-qq')
figname = ['materials8_cemp_fig1_server'];
if print_figs ==1
    cd(figpath)
    export_fig(figname)
    cd(current_dir)
    close
end