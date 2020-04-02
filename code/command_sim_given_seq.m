% command_sim_given_seq.m
% A general way to simulate the model given a sequence of something(s),
% optimizing over input sequences to satisfy residual equations in the
% model.
% 2 April 2020s

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

skip_old_stuff = 1;
%% Initialize

[param, setp, param_names, param_values_str, param_titles] = parameters_next;
sig_r=param.sig_r; sig_i=param.sig_i; sig_u=param.sig_u;
[fyn, fxn, fypn, fxpn] = model_NK(param);
[gx,hx]=gx_hx_alt(fyn,fxn,fypn,fxpn);
[ny, nx] = size(gx);
SIG = eye(nx).*[sig_r, sig_i, sig_u]';
eta = SIG; %just so you know
% Learning ALM matrices
[Aa, Ab, As] = matrices_A_13_true_baseline(param, hx);

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

% Generate one sequence of shocks w/o monpol shock
rng(0)
T = 10
burnin = 0; ne=3;
rng(0)
e = randn(ne,T+burnin);
% turn off monpol shock
e(2,:) = zeros(1,T+burnin);

% An initial simulation
[~, y, ~, ~, ~, ~, ~, ~, ~,~, k] = sim_learnLH(gx,hx,SIG,T+burnin,burnin,e, Aa, Ab, As, param, PLM, gain);
% [~, yi, ~, ~, ~, ~, ~, ~, ~,~, ki] = sim_learnLH_given_i(gx,hx,SIG,T+burnin,burnin,e, Aa, Ab, As, param, PLM, gain, y(3,:));
[~, yg, ~, ~, ~, ~, ~, ~, ~,~, kg] = sim_learnLH_given_seq(gx,hx,SIG,T+burnin,burnin,e, Aa, Ab, As, param, PLM, gain, y(3,:));




%% Optimize over those sequences

%Optimization Parameters
options = optimset('fmincon');
options = optimset(options, 'TolFun', 1e-9, 'display', 'iter');

ub = 40*ones(T,1);
lb = -40*ones(T,1);

% select which of A1-A3 should be a residual equation
s_r = [1;1;1];

%Compute the objective function one time with some values
loss = objective_seq(yg(3,:),param,e,T,burnin,PLM,gain,gx,hx,SIG,Aa,Ab,As,y,s_r)

% return

disp('Begin fmincon... Takes about a minute.')
tic

seq0 = y(3,:);
%Declare a function handle for optimization problem
objh = @(seq) objective_seq(seq,param,e,T,burnin,PLM,gain,gx,hx,SIG,Aa,Ab,As,y,s_r);
[seq_opt, loss_opt] = fmincon(objh, seq0, [],[],[],[],lb,ub,[],options);
% fmincon(FUN,X0,A,B,Aeq,Beq,LB,UB,NONLCON,OPTIONS)
loss
loss_opt
toc

[~, y_opt, ~, ~, ~, ~, ~, ~, ~,~, k_opt] = sim_learnLH_given_seq(gx,hx,SIG,T+burnin,burnin,e, Aa, Ab, As, param, PLM, gain, seq_opt);

y_opt-y

seriesnames = {'\pi', 'x','i'};
figtitle = 'Differences';
create_plot_observables(y_opt-y, seriesnames, figtitle)
figtitle = 'Taylor rule versus optimized';
comparisonnames = {'TR', 'optimized'}; 
create_plot_observables_comparison(y_opt,y, seriesnames, figtitle, comparisonnames)

