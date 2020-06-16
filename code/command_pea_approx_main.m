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
% filename = 'estim_LOMgain_outputs15_Jun_2020';% copo, kmin=0.00001
filename = 'estim_LOMgain_outputs15_Jun_2020_15_45_21'; % copo, kmin=0
load([filename,'.mat'])
xxgrid = estim_outputs{1};
yygrid = estim_outputs{2};
ng     = estim_outputs{3};
k1_opt = estim_outputs{4};
alph_opt = estim_outputs{5};
x = estim_outputs{6};
tol    = estim_outputs{7};
lbname = estim_outputs{8};
ndrop_est  = estim_outputs{9};

alph = alph_opt;

% derivatives of the anchoring function wrt k^{-1}_{t-1} and fe
tic
[g_k,g_fe] = gradient(reshape(k1_opt,ng,ng));
toc
g_pi = g_fe; % input this w/ the k1grid and fegrid so the code doesn't have to take derivatives in every eval
g_pibar = -g_fe;

% fe=0.1;
% kt_1 =1/k1_opt(1);
k1grid = xxgrid(1,:);
fegrid = yygrid(:,1);
% tic
% [k,g_pi,g_pibar] = fk_smooth_approx(alph,x,fe,kt_1, k1grid,fegrid, g_fe)
% toc
% tic
% [k] = fk_smooth_approx(alph,x,fe,kt_1)
% toc



%% Parameters, RE model and Taylor rule
T = 100
N=10; % truly this is 1.

% Some titles for figures
seriesnames = {'\pi', 'x','i'};
invgain = {'Inverse gain'};
residnames = {'IS', 'PC', 'TR'};

% Parameters
ndrop =0; ne=3;
[param, set, param_names, param_values_str, param_titles] = parameters_next;

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
rng(0)
eN = randn(ne,T+ndrop,N);
% zero out the monpol shock
eN(2,:,:) = zeros(T+ndrop,N);
e = squeeze(eN(:,:,1));

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
[x0, y0, k0, phi0, FA0, FB0, FEt_10, diff0] = sim_learnLH_clean_approx(alph,x,param,gx,hx,eta, PLM, gain, T+ndrop,ndrop,e, knowTR,mpshock);
% create_plot_observables(y0,seriesnames, 'Simulation using the Taylor rule', ['implement_anchTC_obs_TR_approx',todays_date], print_figs)
% create_plot_observables(1./k0,invgain, 'Simulation using the Taylor rule', ['implement_anchTC_invgain_TR_approx',todays_date], print_figs)
% return
% %%% when saving to draft or presentations, use these 3 lines below
% cd '/Users/lauragati/Dropbox/BC_Research/next/code'
% create_pretty_subplots(y0,{'$\pi$', '$x$','$i$'}, ['implement_anchTC_obs_TR_approx',todays_date], print_figs)
% create_pretty_plot_x(1:length(k0),1./k0, ['implement_anchTC_invgain_TR_approx',todays_date], print_figs)

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
bet0 = ones(4,3); % 4 states x 3 powers
bet = bet0(:);



% Evaluate residuals once
resids = objective_seq_clean_parametricE_approx(seq0crop,bet,n_inputs,param,gx,hx,eta,PLM,gain,T,ndrop,e, alph, x,  k1grid,fegrid, g_fe, knowTR);
disp('Initial residuals IS, PC, TC and A7')
disp(num2str(resids))

return

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
    objh = @(seq) objective_seq_clean_parametricE_approx(seq,bet,n_inputs,param,gx,hx,eta,PLM,gain,T,ndrop,e, alph,x, k1grid,fegrid, g_fe, knowTR);
    
    tic
    [seq_opt, resids_opt, flag] = fsolve(objh,seq0crop, options);
    toc
    % seq_opt-[seq0(:,2:end-1);FEt_10(1,2:end-1)]
    if flag==1 % If fsolve converged to a root
    % Projection step: Recover v, compute analogues E, update beta
    bet1 = projection_approx(seq_opt,n_inputs,param,gx,hx,eta,PLM,gain,T,ndrop,e,alph,x, k1grid,fegrid, g_fe, knowTR);
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

[~, ysim7, k7, phi7] = sim_learnLH_clean_given_seq3_approx(param,gx,hx,eta,PLM, gain, T,ndrop,e,seq_opt,n_inputs, alph,x, k1grid,fegrid, g_fe);
create_plot_observables(ysim7,seriesnames, 'Simulation using input sequence ', ['implement_anchTC_obs_approx',nowstr], print_figs)
create_plot_observables(1./k7, invgain,'Simulation using input sequence', ['implement_anchTC_invgain_approx', nowstr], print_figs)

if flag==1 % only save output if last fsolve solved.
output = {e, ysim7, k7, phi7, seq_opt};
filename = ['pea_outputs_approx', nowstr, '.mat'];
save(filename, 'output')
disp(['Saving results as ', filename])
end

disp('Done.')