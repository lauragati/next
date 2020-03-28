% command_implement_target_criterion.m
% Numerical implementation of target criterion to get back the optimal
% interest rate path
% developed in materials22.m and replaces that file.
% 26 March 2020

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
%%

[param, set, param_names, param_values_str, param_titles] = parameters_next;

% Structure will be
% 1.) Generate a guess exog interest rate sequence
rng(0)
T = 100; H = 50; burnin = 0; N=1; ne=3;
i_seq0 = gen_AR1(T+H,0.3,1);
% and gen shocks
rng(0)
eN = randn(ne,T+H+burnin,N);

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

% An initial check to see if simulation given i works
sig_r=param.sig_r; sig_i=param.sig_i; sig_u=param.sig_u;
[fyn, fxn, fypn, fxpn] = model_NK(param);
[gx,hx]=gx_hx_alt(fyn,fxn,fypn,fxpn);
[ny, nx] = size(gx);
[Aa, Ab, As] = matrices_A_13_true_baseline(param, hx);
SIG = eye(nx).*[sig_r, sig_i, sig_u]';
eta = SIG; %just so you know
[~, y, ~, ~, ~, ~, ~, ~, ~,~, k] = sim_learnLH(gx,hx,SIG,T+H+burnin,burnin,eN, Aa, Ab, As, param, PLM, gain);
[~, yi, ~, ~, ~, ~, ~, ~, ~,~, k_i] = sim_learnLH_given_i(gx,hx,SIG,T+H+burnin,burnin,eN, Aa, Ab, As, param, PLM, gain, y(3,:));
yi(3,:) - y(3,:)


%%

% 2.) Solve for optimal interest rate path using a target-criterion-based
% loss.
% I'll call this most optimal of plans the Ramsey plan. 

disp('Begin fmincon... Takes about a minute.')
%Optimization Parameters
options = optimset('fmincon');
options = optimset(options, 'TolFun', 1e-9, 'display', 'iter');

ub = 40*ones(T+H,1);
lb = -40*ones(T+H,1);
% %Compute the objective function one time with some values
loss = objective_target_criterion(i_seq0,param,eN,T,N,burnin,PLM,gain);

tic

%Declare a function handle for optimization problem
objh = @(i_seq) objective_target_criterion(i_seq,param,eN,T,N,burnin,PLM,gain);
[i_ramsey, loss_opt] = fmincon(objh, i_seq0, [],[],[],[],lb,ub,[],options);
% fmincon(FUN,X0,A,B,Aeq,Beq,LB,UB,NONLCON,OPTIONS)

toc

% 3.) Simulate model given the Ramsey plan for i -> obtain Ramsey plans for
% x and pi
[pi_ramsey,x_ramsey,k] = fun_sim_anchoring_given_i(param,T+H,N,burnin,eN,PLM,gain,i_ramsey);

% see (plot) deviations between the Ramsey plan and the plans obtained when
% using a Taylor rule
[y_TR,k_TR] = fun_sim_anchoring(param,T+H,N, burnin,eN,PLM,gain);

figure
subplot(2,1,1)
plot(pi_ramsey); hold on
plot(y_TR(1,:))
legend('Ramsey', 'Taylor rule')
title('Inflation')

subplot(2,1,2)
plot(i_ramsey); hold on
plot(y_TR(3,:))
legend('Ramsey', 'Taylor rule')
title('Interest rate')

% 4.) Can even do search over TR parameters to see if they can implement
% the Ramsey plan.