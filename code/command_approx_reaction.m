% Find the interest rate reaction function as an approximation (Taylor- or
% Fourier) of the nonlinear true optimal reaction function. 
% 9 April 2020

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
ne = 3; burnin =0;

T = 40 % length of time
N = 50 % size of cross-section
H = 20  % additional time periods for target criterion

rng(0)
eN = randn(ne,T+H+burnin,N);
% turn off monpol shock
eN(2,:,:) = zeros(T+H+burnin, N);

%% Optimization 1

%Optimization Parameters
options = optimset('fmincon');
options = optimset(options, 'TolFun', 1e-9, 'display', 'iter');


ub = [5,5];
lb = [1.0001,0];
coeffs0 = [1.2,4];

%Compute the objective function one time with some values
tic
loss = objective_approx_reaction(coeffs0,param,eN,T+H,burnin,PLM,gain,gx,hx,SIG,Aa,Ab,As,H,N) % 4 seconds!
toc

disp('Begin fmincon...')
tic

%Declare a function handle for optimization problem
objh = @(coeffs) objective_approx_reaction(coeffs,param,eN,T+H,burnin,PLM,gain,gx,hx,SIG,Aa,Ab,As,H,N);
[coeffs_opt, loss_opt, exitflag,output] = fmincon(objh, coeffs0, [],[],[],[],lb,ub,[],options);
% fmincon(FUN,X0,A,B,Aeq,Beq,LB,UB,NONLCON,OPTIONS)

disp(['psi_pi, psi_x = ', num2str(coeffs_opt), '; for a loss of ', num2str(loss_opt)])

%% Plot loss
psi_x = 0;
psi_pi_vals = linspace(1.001,2,2);
loss = 0;
tic
for psi_pi = psi_pi_vals(1:end)
    loss(end+1) = objective_approx_reaction([psi_pi, psi_x],param,eN,T+H,burnin,PLM,gain,gx,hx,SIG,Aa,Ab,As,H,N);
end
toc