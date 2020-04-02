% command_implement_target_criterion.m
% Numerical implementation of target criterion to get back the optimal
% interest rate path
% developed in materials22.m and replaces that file.
% 26 March 2020
% needs to be reworked once I have found a general way of obtaining
% sequences that work for the model.

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


% Structure will be
% 1.) Generate one sequence of shocks w/o monpol shock
rng(0)
T = 10
H = 1 
burnin = 0; N=1; ne=3;
% sorry, burnin is not possible for this exercise.
rng(0)
e = randn(ne,T+H+burnin);
% turn off monpol shock
e(2,:) = zeros(1,T+H+burnin);

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
[~, y, ~, ~, ~, ~, ~, ~, ~,~, k] = sim_learnLH(gx,hx,SIG,T+H+burnin,burnin,e, Aa, Ab, As, param, PLM, gain);
[~, yi, ~, ~, ~, ~, ~, ~, ~,~, k_i] = sim_learnLH_given_i(gx,hx,SIG,T+H+burnin,burnin,e, Aa, Ab, As, param, PLM, gain, y(3,:));
% yi(1,:) - y(1,:) % these two explode
% yi(2,:) - y(2,:) % these two explode
% yi(3,:) - y(3,:) % these are equal, thats good

% Set the initial Taylor rule implied interest rate path as initial guess
iseq0 = y(3,:);
return
%%

% 2.) Solve for optimal interest rate path using a target-criterion-based
% loss.
% I'll call this most optimal of plans the Ramsey plan. 

%Optimization Parameters
options = optimset('fmincon');
options = optimset(options, 'TolFun', 1e-9, 'display', 'iter');

ub = 40*ones(T+H,1);
lb = -40*ones(T+H,1);
% %Compute the objective function one time with some values
loss = objective_target_criterion(iseq0,param,e,T,burnin,PLM,gain,gx,hx,SIG,Aa,Ab,As)
return


disp('Begin fmincon... Takes about a minute.')
tic

%Declare a function handle for optimization problem
objh = @(iseq) objective_target_criterion(iseq,param,e,T,burnin,PLM,gain,gx,hx,SIG,Aa,Ab,As);
[i_ramsey, loss_opt] = fmincon(objh, iseq0, [],[],[],[],lb,ub,[],options);
% fmincon(FUN,X0,A,B,Aeq,Beq,LB,UB,NONLCON,OPTIONS)
loss
loss_opt
toc



% 3.) Simulate model given the Ramsey plan for i -> obtain Ramsey plans for
% x and pi, and contrast with observables coming from a Taylor rule
[~, y_TR, ~, ~, ~, ~, ~, ~, ~,~, k] = sim_learnLH(gx,hx,SIG,T+H+burnin,burnin,e, Aa, Ab, As, param, PLM, gain);
[~, y_ramsey, ~, ~, ~, ~, ~, ~, ~,~, k_i] = sim_learnLH_given_i(gx,hx,SIG,T+H+burnin,burnin,e, Aa, Ab, As, param, PLM, gain, i_ramsey);


%% Plot deviations between the Ramsey plan and the plans obtained when
% using a Taylor rule

[fs, lw] = plot_configs;

figure
set(gcf,'color','w'); % sets white background color
set(gcf, 'Position', get(0, 'Screensize')); % sets the figure fullscreen
subplot(3,1,1)
plot(y_ramsey(1,:),'linewidth', lw); hold on
plot(y_TR(1,:),'linewidth', lw)
legend('Ramsey', 'Taylor rule', 'location', 'eastoutside')
title('Inflation')
ax = gca; % current axes
ax.FontSize = fs;
grid on
grid minor

subplot(3,1,2)
plot(y_ramsey(2,:),'linewidth', lw); hold on
plot(y_TR(3,:),'linewidth', lw)
legend('Ramsey', 'Taylor rule', 'location', 'eastoutside')
title('Output gap')
ax = gca; % current axes
ax.FontSize = fs;grid on
grid minor

subplot(3,1,3)
plot(i_ramsey,'linewidth', lw); hold on
plot(y_TR(3,:),'linewidth', lw)
legend('Ramsey', 'Taylor rule', 'location', 'eastoutside')
title('Interest rate')
ax = gca; % current axes
ax.FontSize = fs;grid on
grid minor

% 4.) Can even do search over TR parameters to see if they can implement
% the Ramsey plan.