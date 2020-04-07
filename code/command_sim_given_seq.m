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
T = 60 % for T=100 it takes 16 sec; for T=1000 it takes 143 seconds (2.3 minutes)

burnin = 0; ne=3;
rng(0)
e = randn(ne,T+burnin);
% turn off monpol shock
e(2,:) = zeros(1,T+burnin);


% I'm taking T=40 as my benchmark. Since for the anchoring target
% criterion, I cut off the last H=20 periods, if T=60, that means I'm using
% the simple anchoring target criterion
if T==60
    H = 20;
    simple_anchoring_TC = 1;
    % and also replace the last H shocks with the CB's expectation of them
    % conditional on period t information; i.e. zero out innovations
    t=T-H;
    e(:,t:end) = 0;
else
    simple_anchoring_TC = 0;
end

%% Optimize over those sequences

%Optimization Parameters
options = optimset('fmincon');
options = optimset(options, 'TolFun', 1e-9, 'display', 'iter', 'MaxFunEvals', 15000);

%Select exogenous inputs
s_inputs = [1;1;1]; % pi, x, i
i_inputs = find(s_inputs); % index of inputs series in y
n_inputs = sum(s_inputs); % the number of input series

ub = 40*ones(n_inputs,T);
lb = -40*ones(n_inputs,T);

%Compute the objective function one time with some values
loss = objective_seq(rand(n_inputs,T),param,e,T,burnin,PLM,gain,gx,hx,SIG,Aa,Ab,As)

% A simulation using the Taylor rule to intialize
[~, y, ~, ~, ~, ~, ~, ~, ~,~, k] = sim_learnLH(gx,hx,SIG,T+burnin,burnin,e, Aa, Ab, As, param, PLM, gain);
seq0 = y(i_inputs,:);
seq0 = rand(size(seq0));


disp('Begin fmincon... Takes about a minute.')
tic

%Declare a function handle for optimization problem
objh = @(seq) objective_seq(seq,param,e,T,burnin,PLM,gain,gx,hx,SIG,Aa,Ab,As);
[seq_opt, loss_opt] = fmincon(objh, seq0, [],[],[],[],lb,ub,[],options);
% fmincon(FUN,X0,A,B,Aeq,Beq,LB,UB,NONLCON,OPTIONS)
loss
loss_opt
toc

[loss_opt_check, resids_opt] = objective_seq(seq_opt,param,e,T,burnin,PLM,gain,gx,hx,SIG,Aa,Ab,As);

num_res = size(resids_opt,1);
if T==60
    [~, y_opt, ~, ~, ~, ~, ~, ~, ~,~, k_opt] = sim_learnLH_given_seq(gx,hx,SIG,T-H+burnin,burnin,e, Aa, Ab, As, param, PLM, gain, seq_opt);
    [~, y_TR, ~, ~, ~, ~, ~, ~, ~,~, k] = sim_learnLH(gx,hx,SIG,T-H+burnin,burnin,e, Aa, Ab, As, param, PLM, gain);
else
    [~, y_opt, ~, ~, ~, ~, ~, ~, ~,~, k_opt] = sim_learnLH_given_seq(gx,hx,SIG,T+burnin,burnin,e, Aa, Ab, As, param, PLM, gain, seq_opt);
    [~, y_TR, ~, ~, ~, ~, ~, ~, ~,~, k] = sim_learnLH(gx,hx,SIG,T+burnin,burnin,e, Aa, Ab, As, param, PLM, gain);
end
% y_opt-y

%%
seriesnames = {'\pi', 'x','i'};
figtitle = ['Feed in the following exogenous sequences: ', seriesnames{i_inputs}]
comparisonnames = {'TR', 'min residuals'};
% Adapt figname based on whether you imposed TC or not
if num_res==2 % no TC
    disp('No target criterion')
    % Choose whether you initialized at Taylor-rule sequences or at random
    % sequence
    if max(max(abs(seq0-y(i_inputs,:)))) ==0
        figname = strrep([this_code,'_', seriesnames{i_inputs}, '_initialized_atTR'],'\','_');
        disp('Initialized at TR')
    else
        figname = strrep([this_code,'_', seriesnames{i_inputs}, '_initalized_at_rand'],'\','_');
        disp('Initialized at random sequence(s)')
    end
elseif num_res==3 % impose TC
    disp('Imposing target criterion')
    % Choose which target criterion you impose (RE or simple anchoring)
    if simple_anchoring_TC == 0
        disp('RE discretion TC')
        % Choose whether you initialized at Taylor-rule sequences or at random
        % sequence
        if max(max(abs(seq0-y(i_inputs,:)))) ==0
            figname = strrep([this_code,'_', seriesnames{i_inputs}, '_TC_initalized_atTR'],'\','_'); % with TC
            disp('Initialized at TR')
        else
            figname = strrep([this_code,'_', seriesnames{i_inputs}, '_TC_initalized_at_rand'],'\','_'); % with TC
            disp('Initialized at random sequence(s)')
            
        end
    elseif simple_anchoring_TC==1
        disp('simple anchoring TC')
        % Choose whether you initialized at Taylor-rule sequences or at random
        % sequence
        if max(max(abs(seq0-y(i_inputs,:)))) ==0
            figname = strrep([this_code,'_', seriesnames{i_inputs}, '_anch_TC_initalized_atTR'],'\','_'); % with TC
            disp('Initialized at TR')
        else
            figname = strrep([this_code,'_', seriesnames{i_inputs}, '_anch_TC_initalized_at_rand'],'\','_'); % with TC
            disp('Initialized at random sequence(s)')
            
        end
    end
else
    disp('wtf')
end

create_plot_observables_comparison(y_TR,y_opt, seriesnames, figtitle, comparisonnames,figname,print_figs)

disp('Done.')