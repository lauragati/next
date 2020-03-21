% materials22
% 19 March 2020

% Goals:
% 1.) Do initial GMM estimation of the anchoring function in midsimple
% model - think I can take some shortcuts from general codes
% 2.) Do the numerical implementation of the target criterion in midsimple
% model - same shortcuts apply?

clearvars
close all
clc

% Add all the relevant paths and grab the codename
this_code = mfilename;
[current_dir, basepath, BC_researchpath,toolpath,export_figpath,figpath,tablepath,datapath] = add_paths;

% Variable stuff ---
print_figs        = 0;
stop_before_plots = 0;
skip_old_plots    = 0;
output_table = print_figs;

skip_old_stuff = 1;

%% 1.) GMM
% 1.) Simulate anchoring model
[param, set, param_names, param_values_str, param_titles] = parameters_next;
ne = 3;

% Params for the general learning code
constant_only = 1; % learning constant only
constant_only_pi_only = 11; % learning constant only, inflation only
mean_only_PLM = -1;
slope_and_constant = 2;

dgain = 1;  % 1 = decreasing gain, 2 = endogenous gain, 3 = constant gain
again_critCEMP  = 21;
again_critCUSUM = 22;
again_critsmooth = 23;
cgain = 3;

% Model selection
%%%%%%%%%%%%%%%%%%%
PLM = constant_only_pi_only;
gain = again_critsmooth;
%%%%%%%%%%%%%%%%%%%

T = 100 % 400
% Size of cross-section
N = 1 %100, 500
burnin = 0; % 100

% So the simulated data are:
y = fun_sim_anchoring(param,T,N, burnin,ne,PLM,gain);
ny = size(y,1);
% Suppose we have real data (in logs)
y_data = randn(size(y));

% 2.) Filter both using the same filter
% 2.1) HP filter
g = nan(size(y));
c = nan(size(y));
g_data = nan(size(y));
c_data = nan(size(y));
for i=1:ny
[g(i,:),c(i,:)] = HPfilter(y(i,:)');
[g_data(i,:),c_data(i,:)] = HPfilter(y_data(i,:)');
end

% 2.2) Hamilton filter
v = nan(size(y));
v_data = nan(size(y));
for i=1:ny
[v(i,:)] = HPfilter(y(i,:)');
[v_data(i,:)] = HPfilter(y_data(i,:)');
end

% 2.3) BK filter
K=12;
ystar = nan(ny,T-2*K-1);
ystar_data = nan(ny,T-2*K-1);
for i=1:ny
ystar(i,:) = BKfilter(y(i,:)');
ystar_data(i,:) = BKfilter(y(i,:)');
end

% 3.) Get autocovariances
% 4.) Do fmincon