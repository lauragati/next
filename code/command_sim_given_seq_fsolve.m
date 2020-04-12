% command_sim_given_seq_fsolve.m
% As the name says, a version of command_sim_given_seq.m for fsolve.
% 8 April 2020

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

% % Params for the general learning code
% constant_only = 1; % learning constant only
% constant_only_pi_only = 11; % learning constant only, inflation only
% mean_only_PLM = -1;
% slope_and_constant = 2;
% 
% dgain = 1;  % 1 = decreasing gain, 2 = endogenous gain, 3 = constant gain
% again_critCEMP  = 21;
% again_critCUSUM = 22;
% again_critsmooth = 23;
% cgain = 3;

% % Model selection
% %%%%%%%%%%%%%%%%%%%
% PLM = constant_only_pi_only;
% gain = again_critsmooth;
% %%%%%%%%%%%%%%%%%%%


T = 40 

implementTR = 1;
implementRE_TC = 2;
implement_anchTC = 3;

initializeTR = 1;
initialize_rand = 2;

% Specify which optimization to do
%%%%%%%%%%%%%%%%%%%
variant = implementTR;
initialization = initializeTR;
initialization = initialize_rand;
%Select exogenous inputs
s_inputs = [0;0;1]; % pi, x, i
%%%%%%%%%%%%%%%%%%%

H = 0;
ndrop = 0; ne=3;
rng(0)
e = randn(ne,T+ndrop);
% give names for figures
if variant == implementTR
    disp('Variant: Implement Taylor rule')
    variant_title = 'Implement Taylor rule; ';
    variant_name = 'implementTR';
elseif variant == implementRE_TC
    disp('Variant: Implement RE-TC')
    variant_title = 'Implement RE-TC; ';
    variant_name = 'implementRE_TC';
elseif variant == implement_anchTC
    disp('Variant: Implement anchoring-TC')
    variant_title = 'Implement anchoring-TC; ';
    variant_name = 'implement_anchTC';
    H = 20;
    %For the anchoring target criterion, I replace the last H shocks with the CB's expectation of them
    % conditional on period t information; i.e. zero out innovations
    rng(0)
    e = randn(ne,T+H+ndrop);
    e(:,T+1:end) = 0;
end

% turn off monpol shock
e(2,:) = zeros(1,T+H+ndrop);

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

if n_inputs == 3
    disp('Feeding in pi, x, i')
elseif n_inputs==2
    disp('Feeding in x, i')
elseif n_inputs ==1
    disp('Feeding in i')
end


%% Optimize over those sequences

% %Optimization Parameters
options = optimoptions('fsolve', 'TolFun', 1e-9, 'display', 'iter', 'MaxFunEvals', 10000);

% A simulation using the Taylor rule to intialize
[~, y, k_TR, pibar_TR, FA_TR, FB_TR, g_pi_TR, g_pibar_TR, fett_1eve_TR, diff_TR] = sim_learnLH_clean_smooth(param,gx,hx,eta, Aa, Ab, As, T+ndrop,ndrop,e);
seq0 = y(i_inputs,:);
if initialization == initialize_rand
    rng(100)
    seq0 = rand(size(seq0));
end

[xsim, ysim, k, pibar, FA, FB, g_pi, g_pibar, fett_1eve, diff] = sim_learnLH_clean_smooth_given_seq(param,gx,hx,eta,seq0,T+ndrop,ndrop,e);

%Compute the objective function one time with some values
loss = objective_seq_fsolve(seq0,param,gx,hx,eta,T,ndrop,e);
disp('Initial residuals are ')
disp(num2str(loss))

figure
subplot(2,2,1)
plot(loss(1,:)); hold on
plot(loss(2,:))
plot(loss(3,:))
legend('Res TR (A3)', 'Res NKIS (A9)', 'Res NKPC (A10)')
subplot(2,2,2)
plot(diff)
title('Convergence')
subplot(2,2,3)
plot(pibar)
title('pibar')
subplot(2,2,4)
plot(1./k)
title('k^{-1}')


figure
subplot(1,3,1)
plot(pibar,'k'); hold on
plot(pibar_TR, 'b--')
legend('pibar', 'pibar TR')
subplot(1,3,2)
plot(1./k,'k'); hold on
plot(1./k_TR, 'b--')
legend('1/k', '1/k TR')
subplot(1,3,3)
plot(fett_1eve,'k'); hold on
plot(fett_1eve_TR,'b--')
legend('fe_{t|t-1}^{eve}', 'fe_{t|t-1}^{eve} TR')


return

tic
objh = @(seq) objective_seq_fsolve(seq,param,gx,hx,eta,T,ndrop,e);
[seq_opt,FVAL,EXITFLAG, OUTPUT] = fsolve(objh,seq0, options);
toc
% Exitflag values:
% 1  fsolve converged to a root.
% 2  Change in X too small.
% 3  Change in residual norm too small.
% 4  Computed search direction too small.
% 0  Too many function evaluations or iterations.
% -1  Stopped by output/plot function.
% -2  Converged to a point that is not a root.
% -3  Trust region radius too small (Trust-region-dogleg).

message = OUTPUT.message(2:32);
if EXITFLAG == 1
    flag_title = 'FSOLVE: found sol.';
elseif EXITFLAG == 0
    flag_title = 'FSOLVE: stopped prematurely.';
elseif EXITFLAG == -2
    flag_title = 'FSOLVE: no sol.';
else
    disp(EXITFLAG)
    flag_title = ['FSOLVE: ', message ];
end

[~, y_opt] = sim_learnLH_clean_smooth_given_seq(param,gx,hx,eta,seq_opt,T,ndrop,e);
[~, y_TR] = sim_learnLH_clean_smooth(param,gx,hx,eta, Aa, Ab, As, T,ndrop,e);

%% Figures
seriesnames = {'\pi', 'x','i'};
figtitle = [variant_title, 'feeding in: ' seriesnames{i_inputs}, '; \newline ' initialization_title, '; ', flag_title, ...
    '\newline max(max(abs(residual))) = ', num2str(max(max(abs(FVAL))))]
comparisonnames = {'TR', 'min residuals'};
figname = strrep([this_code, '_', variant_name, '_',seriesnames{i_inputs},'_', initialization_name ],'\','_');
% create_plot_observables(y_opt)
create_plot_observables_comparison(y_TR,y_opt, seriesnames, figtitle, comparisonnames,figname,print_figs)

disp('Done.')