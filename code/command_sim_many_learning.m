% command_sim_many_learning
% Simulate many learning models for the LH expectations side-project
% cut together from materials3 and materials10
% 18 Nov 2019

clearvars
close all
clc

% % Add all the relevant paths and grab the codename
this_code = mfilename;
[current_dir, basepath, BC_researchpath,toolpath,export_figpath,figpath,tablepath,datapath] = add_paths;

% Variable stuff ---
print_figs        = 0;
stop_before_plots = 0;
skip_old_plots    = 0;
output_table = print_figs;

skip_old_stuff = 1;

%% Simulation
tic
burnin = 0;

[param, setp] = parameters_next;

bet = param.bet;
sig = param.sig;
alph = param.alph;
kapp = param.kapp;
psi_x = param.psi_x;
psi_pi = param.psi_pi;
w = param.w;
gbar = param.gbar;
thetbar = param.thetbar;
rho_r = param.rho_r;
rho_i = param.rho_i;
rho_u = param.rho_u;
sig_r = param.sig_r;
sig_i = param.sig_i;
sig_u = param.sig_u;
rho = param.rho;
ne = 3;
nx = 4;% now n becomes 4
P = eye(ne).*[rho_r, rho_i, rho_u]';
SIG = eye(nx).*[sig_r, sig_i, sig_u, 0]';


% RE model
[fyn, fxn, fypn, fxpn] = model_NK_intrate_smoothing(param);
[gx,hx]=gx_hx_alt(fyn,fxn,fypn,fxpn);
[ny, nx] = size(gx);
[Ap_RE, As_RE, Aa, Ab, As] = matrices_A_intrate_smoothing(param, hx);


%% Simulate models


% Params for the general learning code
constant_only = 1; % learning constant only
mean_only_PLM = -1;
slope_and_constant = 2;
% lets alternate between these
PLM = slope_and_constant;
dgain = 1;  % 1 = decreasing gain, 2 = endogenous gain, 3 = constant gain
again = 2;
cgain = 3;
critCEMP=1;
critCUSUM=2;

T = 100 % 400
% Size of cross-section
N = 1; %500
eN = randn(ne,T,N); % gen all the N sequences of shocks at once.

for n=1:N
    % Sequence of innovations
    e = [squeeze(eN(:,:,n)); zeros(1,T)]; % adding zero shocks to interest rate lag
    
    % Use Ryan's code to simulate from the RE model
    [~, y_RE] = sim_model(gx,hx,SIG,T,burnin,e);
    
    % LH learning (learning both slope and constant of a vector)
    [~, y_LH] = sim_learnLH(gx,hx,SIG,T,burnin,e, Aa, Ab, As, param, PLM, dgain);
    
    % Euler equation learning (learning both slope and constant). Only this
    % is from materials3.m.
    [~, y_EE] = sim_learn_EE(gx,hx,fxpn,fxn,fypn,fyn,SIG,T,burnin,e,param, PLM, dgain);
end

Y(:,:,1) = y_RE;
Y(:,:,2) = y_LH;
Y(:,:,3) = y_EE;
toc
%% Plots

% only plot last 100 obs
rearranged = permute(Y,[3,2,1]); % swap 3rd and 1st dim for plotting
% series = rearranged;
series = rearranged(:,end-19:end,:); % last 20 obs

figtitle = '';
seriesnames = {'RE', 'LH', 'EE'};
titles_obs = {'Inflation','Output gap','Interest rate'};
figname = [this_code, '_simulatedseries_T', num2str(T)];
subplot_names = titles_obs;
create_subplot(series,subplot_names,figname,print_figs,figtitle,seriesnames)