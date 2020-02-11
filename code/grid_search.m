% grid_search.m
% Search the paramspace for (psi_pi, psi_x) to minimize CB's loss
% Begun in materials14.m
% 26 Jan 2020

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

skip_old_stuff = 1;

%% Parameters
[param, setp, param_names, param_values_str, param_titles] = parameters_next;
ne = 3;

burnin = 0;
T = 400 % 400
% Size of cross-section
N = 100 %500

% Params for the general learning code
constant_only = 1; % learning constant only
mean_only_PLM = -1;
slope_and_constant = 2;

dgain = 1;  % 1 = decreasing gain, 2 = endogenous gain, 3 = constant gain
again_critCEMP  = 21;
again_critCUSUM = 22;
cgain = 3;

% Model selection
%%%%%%%%%%%%%%%%%%%
PLM = constant_only;
gain = again_critCUSUM;
%%%%%%%%%%%%%%%%%%%
[PLM_name, gain_name, gain_title] = give_names(PLM, gain);

disp(['(psi_x, thetbar, thettilde)=   ', num2str([param.psi_x, param.thetbar, param.thettilde])])

% gen all the N sequences of shocks at once.
rng(0)
eN = randn(ne,T,N);

%Optimization Parameters
options = optimset('fmincon');
options = optimset(options, 'TolFun', 1e-9, 'display', 'iter');

%% Fmincon
% takes about 4 min
tic
% varp0 = [1.19,0];
% ub = [1.2,1];
% lb = [1.01,-.01];
varp0 = 1.5;
ub = 2; % 1.5
lb = 1.001;
%Compute the objective function one time with some values
loss = objective_CB(varp0,setp,eN,burnin,PLM,gain);

%Declare a function handle for optimization problem
objh = @(varp) objective_CB(varp,setp,eN,burnin,PLM,gain);
[par_opt, loss_opt] = fmincon(objh, varp0, [],[],[],[],lb,ub,[],options);
% fmincon(FUN,X0,A,B,Aeq,Beq,LB,UB,NONLCON,OPTIONS)
toc

disp(['Learning: psi_pi_opt = ', num2str(par_opt), ', L(psi_pi_opt) = ', num2str(loss_opt)])



%% Fmincon RE
% takes about 4 min
tic

% varp0 = [1.19,0];
% ub = [1.2,1];
% lb = [1.01,-.01];
varp0 = 1.01;
ub = 10; % 1.5
lb = 1.001;
%Compute the objective function one time with some values
loss = objective_CB_RE(varp0,setp,eN,burnin);

%Declare a function handle for optimization problem
objh = @(varp) objective_CB_RE(varp,setp,eN,burnin);
[par_opt_RE, loss_opt_RE] = fmincon(objh, varp0, [],[],[],[],lb,ub,[],options);
% fmincon(FUN,X0,A,B,Aeq,Beq,LB,UB,NONLCON,OPTIONS)
toc

disp(['RE: psi_pi_opt = ', num2str(par_opt_RE), ', L(psi_pi_opt) = ', num2str(loss_opt_RE)])
loss_at_psi_learn = objective_CB_RE(par_opt,setp,eN,burnin);
disp(['RE loss at learning-optimal psi_pi = ', num2str(loss_at_psi_learn)])
