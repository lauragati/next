% materials26
%  prepare a clean code of the implementation of TC for Ryan
% 18 April 2020

clearvars
close all
clc

% Add all the relevant paths and grab the codename
this_code = mfilename;
[current_dir, basepath, BC_researchpath,toolpath,export_figpath,figpath,tablepath,datapath] = add_paths;
todays_date = strrep(datestr(today), '-','_');

% Variable stuff ---
print_figs        = 0;
stop_before_plots = 0;
skip_old_plots    = 0;
output_table = print_figs;

skip = 1;
%% An initial evaluation of loss

T = 100
N = 10
ndrop =0; ne=3;
[param, set, param_names, param_values_str, param_titles] = parameters_next;

sig_r = param.sig_r;
sig_i = param.sig_i;
sig_u = param.sig_u;

% Params for the general learning code
constant_only = 1; % learning constant only (vector learning)
constant_only_pi_only = 11; % learning constant only, inflation only (scalar learning)
slope_and_constant = 2;

dgain = 1;  % 1 = decreasing gain, 2 = endogenous gain, 3 = constant gain
again_critCEMP  = 21;
again_critCUSUM = 22;
again_critsmooth = 23;
cgain = 3;

% Model selection
%%%%%%%%%%%%%%%%%%%
PLM = constant_only;
gain = again_critCUSUM;
%%%%%%%%%%%%%%%%%%%

% display which model you're doing
[PLM_name, gain_name, gain_title] = give_names(PLM, gain);

% RE model
[fyn, fxn, fypn, fxpn] = model_NK(param);
[gx,hx]=gx_hx_alt(fyn,fxn,fypn,fxpn);
[ny, nx] = size(gx);
SIG = eye(nx).*[sig_r, sig_i, sig_u]';
eta = SIG; %just so you know

rng(0)
eN = randn(ne,T+ndrop,N);
e = squeeze(eN(:,:,1));
% zero out the monpol shock
eN(2,:,:) = zeros(T+ndrop,N);

% an initial simulation using the Taylor rule
[xsim, ysim, k, phi_seq, FA, FB, diff] = sim_learnLH_clean(param,gx,hx,eta,PLM, gain, T,ndrop,e);
create_plot_observables(ysim)
create_plot_observables(1./k)

% an initial simulation given an interest rate sequence CONT HERE


% Optimization Parameters
options = optimoptions('fmincon', 'TolFun', 1e-12, 'display', 'iter', 'MaxFunEvals', 10000);

%% fsolve