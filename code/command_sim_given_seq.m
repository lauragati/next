% command_sim_given_seq.m
% A general way to simulate the model given a sequence of something(s),
% optimizing over input sequences to satisfy residual equations in the
% model.
% 2 April 2020

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


T = 40 % for T=100 it takes 16 sec; for T=1000 it takes 143 seconds (2.3 minutes)

% Options:
implementTR = 1;
implementRE_TC = 2;
implement_anchTC = 3;

initializeTR = 1;
initialize_rand = 2;

% Specify which optimization to do
%%%%%%%%%%%%%%%%%%%
variant = implement_anchTC;
% initialization = initializeTR;
initialization = initialize_rand;
%Select exogenous inputs
s_inputs = [1;1;1]; % pi, x, i
%%%%%%%%%%%%%%%%%%%

H = 0;
burnin = 0; ne=3;
rng(0)
e = randn(ne,T+burnin);
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
    e = randn(ne,T+H+burnin);
    e(:,T+1:end) = 0;
end

% turn off monpol shock
e(2,:) = zeros(1,T+H+burnin);

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

%Optimization Parameters
options = optimset('fmincon');
options = optimset(options, 'TolFun', 1e-9, 'display', 'iter', 'MaxFunEvals', 10000);


ub = 40*ones(n_inputs,T+H);
lb = -40*ones(n_inputs,T+H);


% A simulation using the Taylor rule to intialize
[~, y] = sim_learnLH(gx,hx,SIG,T+H+burnin,burnin,e, Aa, Ab, As, param, PLM, gain);
seq0 = y(i_inputs,:);
if initialization == initialize_rand
    seq0 = rand(size(seq0));
end

%Compute the objective function one time with some values
loss = objective_seq(seq0,param,e,T+H,burnin,PLM,gain,gx,hx,SIG,Aa,Ab,As, H, variant)


disp('Begin fmincon... can take up to 30 min for anchoring TC!')
tic

% Use a different objective function for the anchoring TC
if variant == implement_anchTC
loss = objective_seq_anchTC(seq0,param,e,T+H,burnin,PLM,gain,gx,hx,SIG,Aa,Ab,As, H)
objh = @(seq) objective_seq_anchTC(seq,param,e,T+H,burnin,PLM,gain,gx,hx,SIG,Aa,Ab,As, H);
else
objh = @(seq) objective_seq(seq,param,e,T+H,burnin,PLM,gain,gx,hx,SIG,Aa,Ab,As, H, variant);
end

% do it
[seq_opt, loss_opt, exitflag,output] = fmincon(objh, seq0, [],[],[],[],lb,ub,[],options);
% fmincon(FUN,X0,A,B,Aeq,Beq,LB,UB,NONLCON,OPTIONS)
loss
loss_opt
toc

message = output.message(2:10); 
if exitflag == 1
    flag_title = 'FMINCON: found sol.';
elseif exitflag == 0
    flag_title = 'FMINCON: stopped prematurely.';
elseif exitflag == -2
    flag_title = 'FMINCON: no sol.';
else
    disp(['exitflag=',num2str(exitflag)])
    flag_title = ['FMINCON: ', message ];
end


[~, y_opt] = sim_learnLH_given_seq(gx,hx,SIG,T+burnin,burnin,e, Aa, Ab, As, param, PLM, gain, seq_opt, H, variant);
[~, y_TR] = sim_learnLH(gx,hx,SIG,T+burnin,burnin,e, Aa, Ab, As, param, PLM, gain);



%% Figures
seriesnames = {'\pi', 'x','i'};
figtitle = [variant_title, 'feeding in: ' seriesnames{i_inputs}, '; \newline ' initialization_title, '; ', flag_title, ...
    '\newline max(max(abs(residuals))) = ', num2str(loss_opt)]
comparisonnames = {'TR', 'min residuals'};
figname = strrep([this_code, '_', variant_name, '_',seriesnames{i_inputs},'_', initialization_name ],'\','_');

create_plot_observables_comparison(y_TR,y_opt, seriesnames, figtitle, comparisonnames,figname,print_figs)

disp('Done.')

%% Use output of first optimization to reoptimize, now over i only
return
seq0 = seq_opt(end,:); % now initialize i at old opt

s_inputs = [0;0;1]; % now i only
i_inputs = find(s_inputs); % index of inputs series in y
n_inputs = sum(s_inputs); % the number of input series

if n_inputs == 3
    disp('Feeding in pi, x, i')
elseif n_inputs==2
    disp('Feeding in x, i')
elseif n_inputs ==1
    disp('Feeding in i')
end



[seq_opt, loss_opt, exitflag,output] = fmincon(objh, seq0, [],[],[],[],lb,ub,[],options);
% fmincon(FUN,X0,A,B,Aeq,Beq,LB,UB,NONLCON,OPTIONS)
loss
loss_opt
toc

message = output.message(2:10); 
if exitflag == 1
    flag_title = 'FMINCON: found sol.';
elseif exitflag == 0
    flag_title = 'FMINCON: stopped prematurely.';
elseif exitflag == -2
    flag_title = 'FMINCON: no sol.';
else
    disp(['exitflag=',num2str(exitflag)])
    flag_title = ['FMINCON: ', message ];
end


[~, y_opt] = sim_learnLH_given_seq(gx,hx,SIG,T+burnin,burnin,e, Aa, Ab, As, param, PLM, gain, seq_opt, H, variant);
[~, y_TR] = sim_learnLH(gx,hx,SIG,T+burnin,burnin,e, Aa, Ab, As, param, PLM, gain);



%% Figures
seriesnames = {'\pi', 'x','i'};
figtitle = [variant_title, 'feeding in: ' seriesnames{i_inputs}, '; \newline ' initialization_title, '; ', flag_title, ...
    '\newline max(max(abs(residuals))) = ', num2str(loss_opt)]
comparisonnames = {'TR', 'min residuals'};
figname = strrep([this_code, '_', variant_name, '_',seriesnames{i_inputs},'_', initialization_name ],'\','_');

create_plot_observables_comparison(y_TR,y_opt, seriesnames, figtitle, comparisonnames,figname,print_figs)

disp('Done.')