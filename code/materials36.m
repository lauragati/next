% materials36
% Continue to estimate the 1D approximated anchoring function, with
% autocovariogram and different seeds for simulated data and estimation
% 24 June 2020

clearvars
close all
clc

% Add all the relevant paths and grab the codename
this_code = mfilename;
[current_dir, basepath, BC_researchpath,toolpath,export_figpath,figpath,tablepath,datapath, inputsRyan_path] = add_paths;
todays_date = strrep(datestr(today), '-','_');
nowstr = strrep(strrep(strrep(datestr(now), '-','_'), ' ', '_'), ':', '_');

% Variable stuff ---
print_figs        = 0;
stop_before_plots = 0;
skip_old_plots    = 0;
output_table = print_figs;

skip = 1;
[fs, lw] = plot_configs;
datestr(now)

do21 = 0;
do22 = 0;
do23 = 0;
do24 = 0;
do25 = 1;

save_stuff=0;

%% Hopefully general params

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%% Truth (Simulated data) %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

filename = 'acf_sim_univariate_data_21_Jun_2020'; % simulated data, nfe = 6. Note: I'm using the large moments vector.

load([filename, '.mat'])
Om = acf_outputs{1}; % this is the moments vector
Om_boot = acf_outputs{2}; % moments vectors in bootstrapped samples
ny = acf_outputs{3};
p = acf_outputs{4};
K = acf_outputs{5};
filt_data = acf_outputs{6};
lost_periods = acf_outputs{7};
T = length(filt_data);
T = T+lost_periods; % to make up for the loss due to filtering

if contains(filename,'sim')
    alph_true = acf_outputs{8};
    nfe_true  = acf_outputs{9};
end

% Note: 3 moments at lag 0 are repeated. So technically we only have 42
% moments (and they could be correlated further)
% reshape(Om(1:9),3,3)

% return

% weighting matrix for GMM
% note: 2nd argument of var(X,W,DIM) signifies whether to normalize by N-1
% or N. 0 -> N-1, 1 -> N. The third argument, DIM, says along which
% dimension to take the variance.
W = diag(var(Om_boot,0,2));
W1 = W^(-1);

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%% Model %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
[param, setp, param_names, param_values_str, param_titles] = parameters_next;
ne = 3;

sig_r = param.sig_r;
sig_i = param.sig_i;
sig_u = param.sig_u;

% RE model
[fyn, fxn, fypn, fxpn] = model_NK(param);
[gx,hx]=gx_hx_alt(fyn,fxn,fypn,fxpn);
[ny, nx] = size(gx);

% [Aa, Ab, As] = matrices_A_13_true_baseline(param, hx);
SIG = eye(nx).*[sig_r, sig_i, sig_u]';
eta = SIG; %just so you know

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

% display which model you're doing
[PLM_name, gain_name, gain_title] = give_names(PLM, gain);

% Specify info assumption on the Taylor rule and not to include a monpol
% shock
knowTR =1
mpshock=1
%%%%%%%%%%%%%%%%%%%

% Size of cross-section
% we're not doing a whole cross-section here
ndrop = 5 % 0-50

% gen all the N sequences of shocks at once.
rng(1) % rng('default')=rng(0)is the one that was used to generate the true data. 
e = randn(ne,T+ndrop);

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%% Approx %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

nfe=5;
k1min = 0;
k1max=1;
femax = 3.5;
femin = -femax;
fegrid = linspace(femin,femax,nfe);
x = cell(1,1);
x{1} = fegrid;
ng_fine = 100;
fegrid_fine = linspace(femin,femax,ng_fine);