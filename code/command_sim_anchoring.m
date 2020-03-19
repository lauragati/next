% command_sim_anchoring
% Simulate anchoring true-baseline learning models
% adapted from command_IRFs_anchoring.m
% 19 March 2020

clearvars
close all
clc

date_today = strrep(datestr(today, 'yyyy-mm-dd'), '-','_');

% % Add all the relevant paths and grab the codename
this_code = mfilename;
[current_dir, basepath, BC_researchpath,toolpath,export_figpath,figpath,tablepath,datapath] = add_paths;

% Variable stuff ---
print_figs        = 0;
stop_before_plots = 0;
skip_old_plots    = 0;
output_table = print_figs;

plot_IRFs=0;
plot_simulated_sequence = 0;
plot_gains=0;
plot_gain_IRF = 0;
plot_IRFs_anch = 1;
skip_old_stuff = 1;

%% Parameters
tic
[param, set, param_names, param_values_str, param_titles] = parameters_next;


psi_x = param.psi_x;
psi_pi = param.psi_pi;
thetbar =param.thetbar;
thettilde = param.thettilde;

sig_r = param.sig_r;
sig_i = param.sig_i;
sig_u = param.sig_u;
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


%% Model selection and informational assumption

% Model selection
%%%%%%%%%%%%%%%%%%%
PLM = constant_only_pi_only;
gain = again_critsmooth;
%%%%%%%%%%%%%%%%%%%

T = 100 % 400
% Size of cross-section
N = 1 %100, 500
burnin = 0; % 100
h = 10; % h-period IRFs

% RE model
[fyn, fxn, fypn, fxpn] = model_NK(param);
[gx,hx]=gx_hx_alt(fyn,fxn,fypn,fxpn);
[ny, nx] = size(gx);
[Aa, Ab, As] = matrices_A_13_true_baseline(param, hx);
SIG = eye(nx).*[sig_r, sig_i, sig_u]';
eta = SIG; %just so you know

[PLM_name, gain_name, gain_title] = give_names(PLM, gain);


%% Simulate models

% gen all the N sequences of shocks at once.
rng(0)
eN = randn(ne,T+burnin,N);

% Preallocate
y = zeros(ny,T,N);

for n=1:N
    % Sequence of innovations
    e = squeeze(eN(:,:,n));
    
    % dbstop if warning
    % RE
    [x_RE, y_RE] = sim_model(gx,hx,SIG,T,burnin,e);
    % Learning
    [x_LH, y_LH, ~, ~, ~, ~, ~, ~, ~,~, k] = sim_learnLH(gx,hx,SIG,T+burnin,burnin,e, Aa, Ab, As, param, PLM, gain);
    y(:,:,n) = y_LH;
end
k
% y
disp(['(psi_x, psi_pi, thetbar, thettilde)=   ', num2str([psi_x, psi_pi, thetbar, thettilde])])
toc
