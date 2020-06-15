% command_pea
% do parameterized expectations for my model
% takes around 5 min, with adaptive starting point less than 2 min
% migrated from input_sequences_for_Ryan; main_file.m
% 29 May 2020

clearvars
close all
clc

% Add all the relevant paths and grab the codename
this_code = mfilename;
% [current_dir, basepath, BC_researchpath,toolpath,export_figpath,figpath,tablepath,datapath, inputsRyan_path] = add_paths;
todays_date = strrep(datestr(today), '-','_');
nowstr = strrep(strrep(strrep(datestr(now), '-','_'), ' ', '_'), ':', '_');

% Variable stuff ---
print_figs        = 0;
stop_before_plots = 0;
skip_old_plots    = 0;
output_table = print_figs;

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
% Call smat to check what info assumption on the Taylor rule you're using
[s1, s2, s3, s4, s5] = smat(param,hx);
%%%%%%%%%%%%%%%%%%%

% find indeces and number of input sequences
i_inputs = find(s_inputs); % index of inputs series in y
n_inputs = sum(s_inputs); % the number of input series

% display which model you're doing
[PLM_name, gain_name, gain_title] = give_names(PLM, gain);

% an initial simulation using the Taylor rule
[x0, y0, k0, phi0, FA0, FB0, FEt_10, diff0] = sim_learnLH_clean(param,gx,hx,eta,PLM, gain, T,ndrop,e);
create_plot_observables(y0,seriesnames, 'Simulation using the Taylor rule', [this_code, '_implement_anchTC_obs_TR_',PLM_name,'_', todays_date], print_figs)
create_plot_observables(1./k0,invgain, 'Simulation using the Taylor rule', [this_code, '_implement_anchTC_invgain_TR_',PLM_name,'_', todays_date], print_figs)

% % when saving to draft or presentations, use these 3 lines below
% cd '/Users/lauragati/Dropbox/BC_Research/next/code'
% create_pretty_subplots(y0,{'$\pi$', '$x$','$i$'}, [this_code, '_implement_anchTC_obs_TR_',PLM_name,'_', todays_date], print_figs)
% create_pretty_plot_x(1:length(k0),1./k0, [this_code, '_implement_anchTC_invgain_TR_',PLM_name,'_', todays_date], print_figs)




% Note: I'm not inputting anything exogenous for period t=1 b/c that
% just causes errors that by construction fsolve can't close
seq0 = zeros(n_inputs,T);
seq0(:,2:end) = y0(i_inputs,2:end);


%% Parameterized expectations

% %Optimization Parameters
options = optimoptions('fsolve', 'TolFun', 1e-9, 'display', 'iter', 'InitDamping',100000);%, 'MaxFunEvals', 4000);
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
objh = @(seq) objective_seq_clean_parametricE(seq,bet,n_inputs,param,gx,hx,eta,PLM,gain,T,ndrop,e);
% Evaluate residuals once
resids = objective_seq_clean_parametricE(seq0crop,bet,n_inputs,param,gx,hx,eta,PLM,gain,T,ndrop,e);
disp('Initial residuals IS, PC, TC and A7')
disp(num2str(resids))
% return

maxiter=30;
BET = zeros(length(bet),maxiter);
iter=0;
crit=1;
start = now;
datestr(now)
while crit > 1e-6 && iter < maxiter
    iter=iter+1
    BET(:,iter) = bet; % storing betas
    % Now solve model equations given conjectured E
    objh = @(seq) objective_seq_clean_parametricE(seq,bet,n_inputs,param,gx,hx,eta,PLM,gain,T,ndrop,e);
    
    tic
    [seq_opt, resids_opt] = fsolve(objh,seq0crop, options);
    toc
    % seq_opt-[seq0(:,2:end-1);FEt_10(1,2:end-1)]
    % Projection step: Recover v, compute analogues E, update beta
    bet1 = projection(seq_opt,n_inputs,param,gx,hx,eta,PLM,gain,T,ndrop,e);
    crit=max(abs(bet-bet1))
    bet=bet1;
    seq0crop = seq_opt; %<-- this speeds up to 1/3 of the time!!
end
endt = now;
elapsed_seconds = etime(datevec(endt), datevec(start));
disp(['Elapsed: ' num2str(elapsed_seconds), ' sec.'])

[~, ysim7, k7, phi7] = sim_learnLH_clean_given_seq3(param,gx,hx,eta,PLM, gain, T,ndrop,e,seq_opt,n_inputs);
create_plot_observables(ysim7(:,1:end-1),seriesnames, 'Simulation using input sequence ', 'implement_anchTC_obs', print_figs)
create_plot_observables(1./k7(1:end-1), invgain, 'Simulation using input sequence', 'implement_anchTC_invgain', print_figs)

% % when saving to draft or presentations, use these 3 lines below
% cd '/Users/lauragati/Dropbox/BC_Research/next/code'
% create_pretty_subplots(ysim7(:,1:end-1),{'$\pi$', '$x$','$i$'}, [this_code, '_implement_anchTC_obs_anchTC_',PLM_name,'_', todays_date], print_figs)
% create_pretty_plot_x(1:length(k7)-1,1./k7(1:end-1), [this_code, '_implement_anchTC_invgain_anchTC_',PLM_name,'_', todays_date], print_figs)

return
output = {e, ysim7, k7, phi7, seq_opt};
filename = ['pea_outputs_', nowstr, '.mat'];
save(filename, 'output')
disp(['Saving results as ', filename])


disp('Done.')