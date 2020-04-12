% command_numerical_ramsey.m
% Try to see if you can obtain optimal i-sequence by numerically solving
% the Ramsey problem
% 11 April 2020

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
%% Initialize params

[param, setp, param_names, param_values_str, param_titles] = parameters_next;
sig_r=param.sig_r; sig_i=param.sig_i; sig_u=param.sig_u;
[fyn, fxn, fypn, fxpn] = model_NK(param);
[gx,hx]=gx_hx_alt(fyn,fxn,fypn,fxpn);
[ny, nx] = size(gx);
SIG = eye(nx).*[sig_r, sig_i, sig_u]';
eta = SIG; %just so you know
% Learning ALM matrices
[Aa, Ab, As] = matrices_A_13_true_baseline(param, hx);
[Aa_april, Ab_april, As_april] = matrices_A_25_true_baseline(param,hx);
Aa = Aa_april;
Ab = Ab_april;
As = As_april;


T = 40 
N = 10;
ndrop=0; ne=3;
rng(0)
eN = randn(ne,T+ndrop,N);
% turn off monpol shock
eN(2,:,:) = zeros(T+ndrop, N);

initializeTR = 1;
initialize_rand = 2;

% Params for the general learning code
constant_only = 1; % learning constant only (vector - _smooth won't work in this case)
constant_only_pi_only = 11; % learning constant only, inflation only (scalar)
slope_and_constant = 2;

% Model selection
%%%%%%%%%%%%%%%%%%%
PLM = constant_only;
% Specify initialization
%%%%%%%%%%%%%%%%%%%
initialization = initializeTR;
% initialization = initialize_rand;
%Select exogenous inputs
s_inputs = [0;0;1]; % pi, x, i
%%%%%%%%%%%%%%%%%%%

if initialization==initializeTR
    disp('Initialized at Taylor rule sequences')
    initialization_title = 'initialized at Taylor rule sequences';
    initialization_name = 'init_TR';
elseif initialization == initialize_rand
    disp('Initialized at random sequences')
    initialization_title = 'initialized at random sequences';
    initialization_name = 'init_rand';
end

i_inputs = find(s_inputs); % index of inputs series in y
n_inputs = sum(s_inputs); % the number of input series
%% Optimization Parameters
options = optimoptions('fmincon', 'TolFun', 1e-9, 'display', 'iter', 'MaxFunEvals', 10000);

% A simulation using the Taylor rule to intialize
% [~, y] = sim_learnLH_clean_smooth(param,gx,hx,eta, Aa, Ab, As, T+ndrop,ndrop,eN(:,:,1));
[~, y] = sim_learnLH_clean_g(param,gx,hx,eta, Aa, Ab, As,PLM, T+ndrop,ndrop,eN(:,:,1));

seq0 = y(i_inputs,:);
if initialization == initialize_rand
    rng(100)
    seq0 = rand(size(seq0));
end

% % Two silly simulations to see if things actually respond to i, and they
% % do.
% [~, ysim_init] = sim_learnLH_clean_smooth_given_seq(param,gx,hx,eta,seq0,T,ndrop,eN(:,:,1));
% [~, ysim_ones] = sim_learnLH_clean_smooth_given_seq(param,gx,hx,eta,ones(size(seq0)),T,ndrop,eN(:,:,1));

[~, ysim, k, phi_seq, FA, FB, diff] = sim_learnLH_clean_g_given_seq(param,gx,hx,eta,seq0,PLM, T,ndrop,eN(:,:,1));
% create_plot_observables(ysim)
% suppose you input a taylor-rule like thing
seq_thingy = 1.6*y(1,:);
[~, ysim_thingy, ~, ~, ~, ~, diff_thingy] = sim_learnLH_clean_g_given_seq(param,gx,hx,eta,seq_thingy,PLM, T,ndrop,eN(:,:,1));
% it's not stable
create_plot_observables(ysim_thingy)


%Compute the objective function one time with some values
loss = objective_ramsey(seq0,param,gx,hx,eta,ndrop,eN);


return
tic
objh = @(seq) objective_ramsey(seq,param,gx,hx,eta,ndrop,eN);
[seq_opt, loss_opt, exitflag,output] = fmincon(objh, seq0, [],[],[],[],[],[],[],options);
% fmincon(FUN,X0,A,B,Aeq,Beq,LB,UB,NONLCON,OPTIONS)
toc

% found local min after 98 sec, 110 iter.

% just choose one sequence of shocks for now, but should take
% cross-sectional average
% [~, y_opt] = sim_learnLH_clean_smooth_given_seq(param,gx,hx,eta,seq_opt,T,ndrop,eN(:,:,1));
% [~, y_TR] = sim_learnLH_clean_smooth(param,gx,hx,eta, Aa, Ab, As, T,ndrop,eN(:,:,1));
[~, y_TR] = sim_learnLH_clean_g(param,gx,hx,eta, Aa, Ab, As,PLM, T+ndrop,ndrop,eN(:,:,1));
[~, y_opt,k_opt, phi_seq_opt, FA_opt, FB_opt, diff_opt] = sim_learnLH_clean_g_given_seq(param,gx,hx,eta,seq_opt,PLM, T,ndrop,eN(:,:,1));


%% Figures
seriesnames = {'\pi', 'x','i'};
figtitle = ['Ramsey policy, ', 'feeding in: ' seriesnames{i_inputs}, '; ', initialization_title]
comparisonnames = {'TR', 'Ramsey policy'};
figname = strrep([this_code, '_',seriesnames{i_inputs},'_', initialization_name ],'\','_');
% create_plot_observables(y_opt)
create_plot_observables_comparison(y_TR,y_opt, seriesnames, figtitle, comparisonnames,figname,print_figs)

disp('Done.')
