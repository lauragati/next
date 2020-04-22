% Main file for inputting exogenous sequence(s) to make all the equations
% of the model hold
% 20 April 2020

clearvars
close all
clc

%% 1.) Parameters - don't need to edit anything here
T = 100
N = 10

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
e = squeeze(eN(:,:,1));
% zero out the monpol shock
eN(2,:,:) = zeros(T+ndrop,N);

% Params for the general learning code
constant_only = 1; % learning constant only (vector learning)
constant_only_pi_only = 11; % learning constant only, inflation only (scalar learning)
slope_and_constant = 2; % (vector learning)

dgain = 1;  % 1 = decreasing gain, 2 = endogenous gain, 3 = constant gain
egain_critCEMP  = 21;
egain_critCUSUM = 22;
egain_critsmooth = 23;
cgain = 3;

init_TR = 1;
init_rand =2;

%%%%%%%%%%%%%%%%%%%
%% 2.) Model selection - this is what you can play around with
%%%%%%%%%%%%%%%%%%%
PLM = constant_only;
gain = egain_critCUSUM;
%%%%%%%%%%%%%%%%%%%
% Initialize input sequence(s)
initialization = init_TR;
% initialization = init_rand;
%%%%%%%%%%%%%%%%%%%
%Select exogenous inputs
s_inputs = [0;0;1]; % pi, x, i
%%%%%%%%%%%%%%%%%%%
% Call smat to check what info assumption on the Taylor rule you're using
[s1, s2, s3, s4, s5] = smat(param,hx);
%%%%%%%%%%%%%%%%%%%

%% 3.) An initial evaluation of loss

% find indeces and number of input sequences
i_inputs = find(s_inputs); % index of inputs series in y
n_inputs = sum(s_inputs); % the number of input series

% display which model you're doing
[PLM_name, gain_name, gain_title] = give_names(PLM, gain);
% display initialization settings
[init_name, init_title] = disp_init_seq(initialization, n_inputs);


% an initial simulation using the Taylor rule
[x0, y0, k0, phi0, FA0, FB0, diff0] = sim_learnLH_clean(param,gx,hx,eta,PLM, gain, T,ndrop,e);
create_plot_observables(y0,seriesnames, 'Simulation using the Taylor rule')
create_plot_observables(1./k0,invgain, 'Simulation using the Taylor rule')


% Note: I'm not inputting anything exogenous for period t=1 b/c that
% just causes errors that by construction fsolve can't close
seq0 = zeros(n_inputs,T);
seq0(:,2:end) = y0(i_inputs,2:end);
if initialization == init_rand
    rng(100)
    seq0(:,2:end) = rand(n_inputs,T-1);
end

% an initial simulation given exogenous input sequence(s) 
[x1, y1, k1, phi1, FA1, FB1, diff1] = sim_learnLH_clean_given_seq(param,gx,hx,eta,PLM, gain, T,ndrop,e,seq0);
create_plot_observables(y1,seriesnames, 'Simulation using input sequence ')
create_plot_observables(1./k1, invgain,'Simulation using input sequence')

% An initial evaluation of objective function
resids = objective_seq_clean(seq0,param,gx,hx,eta,PLM,gain,T,ndrop,e);
disp('Initial residuals are NKIS, NKPC, TR:')
disp(num2str(resids))

create_plot_observables(resids,residnames, 'Errors - Simulation using input sequence ')