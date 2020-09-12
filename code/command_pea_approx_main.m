% command_pea_approx_main
% do parameterized expectations for my model, with outputs of estimation of
% approximated LOM gain / anchoring function
% takes around 5 min
% 16 June 2020

clearvars
close all
clc

% Add all the relevant paths and grab the codename
this_code = mfilename;
[current_dir, basepath,BC_researchpath,toolpath,export_figpath,figpath,tablepath,datapath,tryouts_path] = add_paths;  
todays_date = strrep(datestr(today), '-','_');
nowstr = strrep(strrep(strrep(datestr(now), '-','_'), ' ', '_'), ':', '_');

% Variable stuff ---
print_figs        = 0;
stop_before_plots = 0;
skip_old_plots    = 0;
output_table = print_figs;

% Load estimation outputs

% % filename = 'best_n100_29_Jun_2020'; % materials35 candidate
% % load the saved stuff
% load([filename,'.mat'])
% alph_best = output{1};
% resnorm = output{2};
% alph_opt = alph_best(:,1);
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
% ng_fine = 100;
% fegrid_fine = linspace(femin,femax,ng_fine);

filename = 'estim_LOMgain_outputs_univariate_coax11_Sep_2020_15_46_40'; % materials44 candidate
% load the saved stuff
load([filename,'.mat'])
% Structure of saved file:
% estim_configs={nfe,gridspacing,femax,femin,ub,lb,Wprior,Wdiffs2,Wmid,Wmean,T,ndrop,N,eN, rngsetting};
% learn_configs = {param,PLM_name, gain_name, knowTR, mpshock};
% estim_outputs = {fegrid_fine, ng_fine, alph_opt, alph_k, ALPH0, x, estim_configs, learn_configs};
fegrid_fine = estim_outputs{1};
ng_fine     = estim_outputs{2};
alph_opt      = estim_outputs{3};
alph_k        = estim_outputs{4};
ALPH0         = estim_outputs{5};
x             = estim_outputs{6};
estim_configs = estim_outputs{7};
learn_configs = estim_outputs{8};
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

fegrid = x{1};


% fegrid_uneven = x{1};
% fegrid = fegrid_uneven;
% If you wanna use the uniform grid, then uncomment the following 3 lines:
% fegrid = linspace(femin,femax,nfe);
% x = cell(1,1);
% x{1} = fegrid;

% evaluate gradients of estimated LOM gain beforehand
k1_opt = ndim_simplex_eval(x,fegrid_fine,alph_opt);
g_fe = gradient(k1_opt);


alph = alph_opt;


% return

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% %% 23 August 2020: just try the calibrated values in command_simgas.m (Materials 42)
% 
% alph = [1.0000    0.5000         0    0.5000    1.0000]'
% fegrid = [-4,-3,0,3,4]
% x{1} = fegrid;
% 
% [param, set, param_names, param_values_str, param_titles] = parameters_next;
% 
% sig_r = 0.01;
% sig_i = 2;
% sig_u = 0.5;
% 
% eta = eye(3).*[sig_r, sig_i, sig_u]'

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% %% 27 August 2020: calibration C (Materials 43)
% 
% alph = [0.8    0.4         0    0.4    0.8]'
% fegrid = [-4,-3,0,3,4]
% x{1} = fegrid;
% 
% [param, setp, param_names, param_values_str, param_titles] = parameters_next;

%% Parameters, RE model and Taylor rule
T = 100
N=1; % truly this is 1.

% Some titles for figures
seriesnames = {'\pi', 'x','i'};
invgain = {'Inverse gain'};
residnames = {'IS', 'PC', 'TR'};

% Parameters
ndrop =0; ne=3;
% [param, set, param_names, param_values_str, param_titles] = parameters_next; % now using the ones from the estimation

sig_r = param.sig_r;
sig_i = param.sig_i;
sig_u = param.sig_u;

% RE model
[fyn, fxn, fypn, fxpn] = model_NK(param);
[gx,hx]=gx_hx_alt(fyn,fxn,fypn,fxpn);
[ny, nx] = size(gx);
SIG = eye(nx).*[sig_r, sig_i, sig_u]';
eta = SIG; %just so you know

% Generate innovations
rng(2) % rng(2) default
e = randn(ne,T+ndrop);
% zero out the monpol shock
e(2,:) = zeros(1,T+ndrop);
v = zeros(ne+1,T+ndrop);

% Params for the general learning code
constant_only = 1; % learning constant only (vector learning)
constant_only_pi_only = 11; % learning constant only, inflation only (scalar learning)
slope_and_constant = 2; % (vector learning)

dgain = 1;  % 1 = decreasing gain, 2 = endogenous gain, 3 = constant gain
egain_critCEMP  = 21;
egain_critCUSUM = 22;
egain_critsmooth = 23;
cgain = 3;

%%%%%%%%%%%%%%%%%%%
% Model selection
%%%%%%%%%%%%%%%%%%%
PLM = constant_only_pi_only;
gain = egain_critsmooth;
%%%%%%%%%%%%%%%%%%%
%Select exogenous inputs
s_inputs = [1;1;1]; % pi, x, i
%%%%%%%%%%%%%%%%%%%
% Specify info assumption on the Taylor rule and not to include a monpol
% shock
knowTR =0
mpshock=0
%%%%%%%%%%%%%%%%%%%


% find indeces and number of input sequences
i_inputs = find(s_inputs); % index of inputs series in y
n_inputs = sum(s_inputs); % the number of input series

% display which model you're doing
[PLM_name, gain_name, gain_title] = give_names(PLM, gain);

% an initial simulation using the Taylor rule
disp('Initial simul using TR - this should explode b/c agents don''t know the TR')
[x0, y0, k0, phi0, FA0, FB0, FEt_10, diff0] = sim_learnLH_clean_approx_univariate(alph,x,param,gx,hx,eta, PLM, gain, T+ndrop,ndrop,e,v, knowTR,mpshock);

% % Annualize pi and i
% y0([1,3],:) = ((y0([1,3],:)/100+1).^4 -1)*100;
% 
% % create_plot_observables(y0,seriesnames, 'Simulation using the Taylor rule', ['implement_anchTC_obs_TR_approx',todays_date], print_figs)
% % create_plot_observables(1./k0,invgain, 'Simulation using the Taylor rule', ['implement_anchTC_invgain_TR_approx',todays_date], print_figs)
% % return
% % %%% when saving to draft or presentations, use these 5 lines below
% xlab = 'Quarters';
% ylab = 'dont display';
% xlplus = [45,0.01];
% ylplus = [-400,0];
% create_pretty_subplots(y0,{'$\pi$', '$x$','$i$'},xlab,ylab,xlplus, ylplus, ['implement_anchTC_obs_TR_approx',todays_date], print_figs)
% % create_pretty_plot_x(1:length(k0),1./k0, ['implement_anchTC_invgain_TR_approx',todays_date], print_figs)
% 
% return


% Note: I'm not inputting anything exogenous for period t=1 b/c that
% just causes errors that by construction fsolve can't close
seq0 = zeros(n_inputs,T);
seq0(:,2:end) = y0(i_inputs,2:end);



%% Parameterized expectations

% %Optimization Parameters
options = optimoptions('fsolve', 'TolFun', 1e-9, 'display', 'iter');%, 'MaxFunEvals', 4000);
options.UseParallel=true;
% initDamping = initial value of Levenberg-Marquardt lambda.

% 7.) Parameterized expectations approach
% initialize at Taylor rule
seq0crop = [seq0(:,2:end-1);FEt_10(1,2:end-1)]; % input jumps and Fe(pi)
% Or detach from Taylor rule
seq0crop = rand(size(seq0crop));
% initialize beta-coefficients
bet0 = zeros(4,3); % 4 states x 3 powers % INITIALIZE AT ZEROS INSTEAD OF ONES AND IT FINDS SOL
bet = bet0(:);

% return

% Evaluate residuals once
resids = objective_seq_clean_parametricE_approx(seq0crop,bet,n_inputs,param,gx,hx,eta,PLM,gain,T,ndrop,e, alph, x,fegrid, g_fe, knowTR);
disp('Initial residuals IS, PC, TC and A7')
disp(num2str(resids))

% return

dbstop if warning

maxiter=12;
BET = zeros(length(bet),maxiter);
iter=0;
crit=1;
start = now;
datestr(now)
while crit > 1e-6 && iter < maxiter
    iter=iter+1
    BET(:,iter) = bet; % storing betas
    % Now solve model equations given conjectured E 
    objh = @(seq) objective_seq_clean_parametricE_approx(seq,bet,n_inputs,param,gx,hx,eta,PLM,gain,T,ndrop,e, alph,x,fegrid, g_fe, knowTR);
    tic
    [seq_opt, resids_opt, flag] = fsolve(objh,seq0crop, options);
    toc
    % seq_opt-[seq0(:,2:end-1);FEt_10(1,2:end-1)]
    if flag==1 % If fsolve converged to a root
    % Projection step: Recover v, compute analogues E, update beta
    bet1 = projection_approx(seq_opt,n_inputs,param,gx,hx,eta,PLM,gain,T,ndrop,e,alph,x,fegrid, g_fe, knowTR);
    crit=max(abs(bet-bet1))
    bet=bet1;
    seq0crop = seq_opt; % try to accelerate 
    else
        disp('No sol, stopping PEA.')
        return
    end
end
endt = now;
elapsed_seconds = etime(datevec(endt), datevec(start));
disp(['Elapsed: ' num2str(elapsed_seconds), ' sec.'])

[~, ysim7, k7, phi7] = sim_learnLH_clean_given_seq3_approx(param,gx,hx,eta,PLM, gain, T,ndrop,e,seq_opt,n_inputs, alph,x,fegrid, g_fe, knowTR);
create_plot_observables(ysim7,seriesnames, 'Simulation using input sequence ', ['implement_anchTC_obs_approx',nowstr], print_figs)
create_plot_observables(1./k7, invgain,'Simulation using input sequence', ['implement_anchTC_invgain_approx', nowstr], print_figs)

if flag==1 % only save output if last fsolve solved.
output = {e, ysim7, k7, phi7, seq_opt};
filename = ['pea_outputs_approx', nowstr, '.mat'];
save(filename, 'output')
disp(['Saving results as ', filename])
end

disp('Done.')