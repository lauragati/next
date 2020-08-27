% grid_search_approx.m
% a version of grid_search.m that searches for optimal TR coefficients for
% the case of the univariate approximated LOM gain
% Search the paramspace for (psi_pi, psi_x) to minimize CB's loss
% Begun in materials35.m
% 30 June 2020

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

% Load estimated LOM gain coefficents
% filename = 'best_n100_29_Jun_2020'; % materials35 candidate
% load([filename,'.mat'])
% alph_best = output{1};
% resnorm = output{2};
% alph = alph_best(:,1);
% % grab the rest from materials35, part 2.5
% nfe=5;
% k1min = 0;
% k1max= 1;
% femax = 3.5;
% femin = -femax;
% % and from materials35, intro
% fegrid = linspace(femin,femax,nfe);
% x = cell(1,1);
% x{1} = fegrid;

filename = 'estim_LOMgain_outputs_univariate16_Jul_2020_15_25_10'; % materials37 candidate
load([filename,'.mat'])
% Structure of saved file:
% estim_configs={nfe,gridspacing,femax,femin,ub,lb,Wprior,Wdiffs2,Wmid,Wmean,T,ndrop,N,eN, rngsetting};
% learn_configs = {param, PLM_name, gain_name, knowTR, mpshock};
% estim_outputs = {fegrid_fine, ng_fine, k1_opt, alph_opt_mean, x, estim_configs, learn_configs};
fegrid_fine = estim_outputs{1};
ng_fine     = estim_outputs{2};
k1_opt      = estim_outputs{3};
alph_opt_mean = estim_outputs{4};
x             = estim_outputs{5};
estim_configs = estim_outputs{6};
learn_configs = estim_outputs{7};
nfe            = estim_configs{1};
gridspacing    = estim_configs{2};
femax          = estim_configs{3};
femin          = estim_configs{4};
ub             = estim_configs{5};
lb             = estim_configs{6};
Wprior         = estim_configs{7};
Wdiffs2        = estim_configs{8};
Wmid           = estim_configs{9};
Wmean          = estim_configs{10};
T_est          = estim_configs{11};
ndrop_est      = estim_configs{12};
N_est          = estim_configs{13};
eN_est         = estim_configs{14};
rngsetting_est = estim_configs{15};
param       = learn_configs{1};
PLM_name    = learn_configs{2};
gain_name   = learn_configs{3};
knowTR_est  = learn_configs{4};
mpshock_est = learn_configs{5};

% return
fegrid_uneven = x{1};
fegrid = fegrid_uneven;
% % If you wanna use the uniform grid, then uncomment the following 3 lines:
% fegrid = linspace(femin,femax,nfe);
% x = cell(1,1);
% x{1} = fegrid;

alph = alph_opt_mean;




[param, setp, param_names, param_values_str, param_titles] = parameters_next; % use ones from estimation!
ne = 3;

burnin = 0;
T = 400 % 400
% Size of cross-section
N = 100 %500

% Params for the general learning code
constant_only = 1; % learning constant only
constant_only_pi_only = 11; % learning constant only, inflation only (scalar learning)
slope_and_constant = 2;

dgain = 1;  % 1 = decreasing gain, 2 = endogenous gain, 3 = constant gain
again_critCEMP  = 21;
again_critCUSUM = 22;
again_smooth = 23;

cgain = 3;

% Model selection
%%%%%%%%%%%%%%%%%%%
PLM = constant_only_pi_only;
gain = again_smooth;
%%%%%%%%%%%%%%%%%%%
[PLM_name, gain_name, gain_title] = give_names(PLM, gain);

disp(['(psi_x, lam_x, lam_i)=   ', num2str([param.psi_x, param.lamx, param.lami])])

knowTR=1

% gen all the N sequences of shocks at once.
rng(0)
eN = randn(ne,T,N);

%Optimization Parameters
options = optimset('fmincon');
options = optimset(options, 'TolFun', 1e-9, 'display', 'iter');
options.UseParallel=true;

% %% 23 August 2020: just try the calibrated values in command_simgas.m (Materials 42)
% 
% alph = [1.0000    0.5000         0    0.5000    1.0000]'
% fegrid = [-4,-3,0,3,4]
% x{1} = fegrid;
% 
% [param, setp, param_names, param_values_str, param_titles] = parameters_next;

%% 27 August 2020: calibration C (Materials 43)

alph = [0.8    0.4         0    0.4    0.8]'
fegrid = [-4,-3,0,3,4]
x{1} = fegrid;

[param, setp, param_names, param_values_str, param_titles] = parameters_next;

setp.lamx  = 0.05;
setp.lami  = 0;


%% Fmincon
% dbstop if caught error
% takes about 3 min
tic
% varp0 = [1.19,0];
% ub = [1.2,1];
% lb = [1.01,-.01];
varp0 = 1.5;
ub = 2; % 1.5
lb = 1;
%Compute the objective function one time with some values
loss = objective_CB_approx(varp0,setp,eN,burnin,PLM,gain, alph,x, knowTR);

% return
%Declare a function handle for optimization problem
objh = @(varp) objective_CB_approx(varp,setp,eN,burnin,PLM,gain, alph,x, knowTR);
[par_opt, loss_opt] = fmincon(objh, varp0, [],[],[],[],lb,ub,[],options);
% fmincon(FUN,X0,A,B,Aeq,Beq,LB,UB,NONLCON,OPTIONS)
toc

disp(['Learning: psi_pi_opt = ', num2str(par_opt)])
disp('%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%')



%% Fmincon RE
% takes about 30 sec
tic

% varp0 = [1.19,0];
% ub = [1.2,1];
% lb = [1.01,-.01];
% varp0 = 1.5;
% do separate bounds for RE
ub = 10; % 1.5
lb = 1;
%Compute the objective function one time with some values
loss = objective_CB_RE(varp0,setp,eN,burnin);

%Declare a function handle for optimization problem
objh = @(varp) objective_CB_RE(varp,setp,eN,burnin);
[par_opt_RE, loss_opt_RE] = fmincon(objh, varp0, [],[],[],[],lb,ub,[],options);
% fmincon(FUN,X0,A,B,Aeq,Beq,LB,UB,NONLCON,OPTIONS)
toc

disp(['RE: psi_pi_opt = ', num2str(par_opt_RE)])


return
%% Loss at psi_pi = 1.5, RE- and Anchoring-optimal psi_pi for the baseline parameters of lamx = 0.05, lami = 0

% psi_pi=1.5
RE_loss_at15   = objective_CB_RE(1.5,setp,eN,burnin)
Anch_loss_at15 = objective_CB_approx(1.5,setp,eN,burnin,PLM,gain, alph,x, knowTR)

% psi_pi = RE-optimal
RE_loss_at_REopt   = objective_CB_RE(par_opt_RE,setp,eN,burnin)
Anch_loss_at_REopt = objective_CB_approx(par_opt_RE,setp,eN,burnin,PLM,gain, alph,x, knowTR)

% psi_pi = anch-optimal
RE_loss_at_Anchopt   = objective_CB_RE(par_opt,setp,eN,burnin)
Anch_loss_at_Anchopt = objective_CB_approx(par_opt,setp,eN,burnin,PLM,gain, alph,x, knowTR)


%% Loss under optimal policy:

% pea_output_name = 'pea_outputs_approx23_Aug_2020_14_38_14'; % 23 August 2020 calibration 
pea_output_name = 'pea_outputs_approx27_Aug_2020_14_45_53';  % Calibration C of Materials 43; rng(2) default

load([pea_output_name, '.mat'])
ysim7 = output{2};

% Really we'd need a cross-section of PEA results to invoke the loss
% function:
% EL = loss_CB(setp,ysim7)

% I'll do this manually here for the single history we have available
X2 = ysim7(2,:).^2;
PI2 = ysim7(1,:).^2;
I2  = ysim7(3,:).^2;
bet = setp.bet;
bett = bet.^(1:100);
period_loss = bett .* (PI2 + setp.lamx*X2 + setp.lami*I2);
loss = sum(period_loss)

% Take care of that in PEA we have 100 periods, while here we considered 400
bett2 = bet.^(1:T);
period_loss = bett2 .* repmat(PI2 + setp.lamx*X2 + setp.lami*I2,1,4);
loss = sum(period_loss);
estimate_loss_optimal_policy = loss
