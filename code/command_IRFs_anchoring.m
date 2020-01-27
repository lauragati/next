% command_IRFs_anchoring
% Do IRFs for nachoring true-baseline learning models
% adapted from command_IRFs_many_learning.m
% 25 Jan 2020

clearvars
close all
clc

date_today = strrep(datestr(today, 'yyyy-mm-dd'), '-','_');

% % Add all the relevant paths and grab the codename
this_code = mfilename;
[current_dir, basepath, BC_researchpath,toolpath,export_figpath,figpath,tablepath,datapath] = add_paths;

% Variable stuff ---
print_figs        = 0;
stop_before_plots = 0;
skip_old_plots    = 0;
output_table = print_figs;

plot_IRFs=0;
plot_simulated_sequence = 0;
plot_gains=1;
skip_old_stuff = 1;

%% Parameters
tic
burnin = 0;

[param, set, param_names, param_values_str, param_titles] = parameters_next;

% bet = param.bet;
sig = param.sig;
% alph = param.alph;
% kapp = param.kapp;
psi_x = param.psi_x;
psi_pi = param.psi_pi;
% w = param.w;
% gbar = param.gbar;
% thetbar = param.thetbar;
% rho_r = param.rho_r;
% rho_i = param.rho_i;
% rho_u = param.rho_u;
sig_r = param.sig_r;
sig_i = param.sig_i;
sig_u = param.sig_u;
rho = param.rho;
% p11 = param.p11;
% p22 = param.p22;
ne = 3;

% Params for the general learning code
constant_only = 1; % learning constant only
mean_only_PLM = -1;
slope_and_constant = 2;

dgain = 1;  % 1 = decreasing gain, 2 = endogenous gain, 3 = constant gain
again_critCEMP  = 21;
again_critCUSUM = 22;
cgain = 3;


%% Model selection and informational assumption

% Model selection
%%%%%%%%%%%%%%%%%%%
PLM = constant_only;
gain = again_critCUSUM;
%%%%%%%%%%%%%%%%%%%

T = 400 % 400
% Size of cross-section
N = 100 %500
dt_vals = 25; % time of imposing innovation
h = 10; % h-period IRFs

% RE model
[fyn, fxn, fypn, fxpn] = model_NK(param);
[gx,hx]=gx_hx_alt(fyn,fxn,fypn,fxpn);
[ny, nx] = size(gx);
[Aa, Ab, As] = matrices_A_13_true_baseline(param, hx);
SIG = eye(nx).*[sig_r, sig_i, sig_u]';
eta = SIG; %just so you know

[PLM_name, gain_name, gain_title] = give_names(PLM, gain);


%% Simulate models

% gen all the N sequences of shocks at once.
rng(0)
eN = randn(ne,T,N);

% Preallocate
nd = size(dt_vals,2);
d = 1; % the innovation, delta
GIR_Y_EE = zeros(ny,h,N,nd);
GIR_Y_LH = zeros(ny,h,N,nd);
k = zeros(T,N);
% warning off
for s=2  %2->zoom in on monetary policy shock
    x0 = zeros(1,nx);
    x0(s) = d;
    
    for n=1:N
        % Sequence of innovations
        e = squeeze(eN(:,:,n));
        
        % Unshocked
        dbstop if warning
        % RE
        [x_RE, y_RE] = sim_model(gx,hx,SIG,T,burnin,e);
        % Learning
        [x_LH, y_LH, ~, ~, ~, ~, ~, ~, ~,~, k(:,n)] = sim_learnLH(gx,hx,SIG,T,burnin,e, Aa, Ab, As, param, PLM, gain);
        
        % Shocked
        % RE
        % make RE shock the same scale as learning:
        x0RE = (SIG*x0')';
        [IR, iry, irx]=ir(gx,hx,x0RE,h);
        iry = iry';
        RE_fcsts = gx*hx*irx';
        
        % Learning
        for t=1:nd
            dt = dt_vals(t);
            
            % Shocked
            [~, ys_LH] = sim_learnLH(gx,hx,SIG,T,burnin,e, Aa, Ab, As, param, PLM, gain, dt, x0);
            % Construct GIRs
            GIR_Y_LH(:,:,n,t) = ys_LH(:,dt:dt+h-1) - y_LH(:,dt:dt+h-1);
        end
        
    end
end
% warning on
% Construct RIRs by simple method: means (Option 1)
RIR_Y_LH = squeeze(mean(GIR_Y_LH,3));

disp(['(psi_x, psi_pi, rho, sig)=   ', num2str([psi_x, psi_pi, rho, sig])])
toc

if stop_before_plots==1
    return
end


%% Plots

shocknames = {'natrate', 'monpol','costpush'};
titles_obs = {'Inflation','Output gap','Int. rate'};
titles_fcsts = {'E^m_t(\pi_{t+1})', 'E^e_t(\pi_{t+1})'};
titles_LH = {'fa', 'fb'};
titles_FEs = {'FE^m_t(\pi_{t+1})', 'FE^e_t(\pi_{t+1})'};

% create indexes for the positions of the parameters
fn = fieldnames(param);
make_index(fn)
interesting_param_names = param_names([psi_pi_idx, psi_x_idx, gbar_idx, thetbar_idx, thettilde_idx, kap_idx, alph_cb_idx]);
interesting_param_vals = param_values_str([psi_pi_idx, psi_x_idx, gbar_idx, thetbar_idx, thettilde_idx, kap_idx, alph_cb_idx]);
param_names_vals = cell(size(interesting_param_vals));
relevant_params = 'params';
for i=1:size(param_names_vals,2)
    param_names_vals{i} = [interesting_param_names{i},'_',interesting_param_vals{i}];
    relevant_params = [relevant_params, '_', param_names_vals{i}];
end

if plot_IRFs==1
    for t=1:nd % for the two diff times of imposing the shock
        dt = dt_vals(t);
        
        clear series
        % 1) IRF: OBSERVABLES LH against RE
        series(1,:,:) = RIR_Y_LH(:,:,t)';
        series(2,:,:) = iry';
        figname = [this_code, '_', 'RIR_LH_' shocknames{s}, '_', gain_name, '_', PLM_name , '_', date_today];
        subplot_names = titles_obs;
        legendnames = {'Learning', 'RE'};
        figtitle = [gain_title, '; shock imposed at t=', num2str(dt_vals(t))];
        create_subplot(series,subplot_names,figname,print_figs, figtitle, legendnames)
    end
end

if plot_simulated_sequence==1
    clear series
    T_t = 100;
    % 2) SIMULATED HISTORY: OBSERVABLES LH against RE
    series(1,:,:) = y_LH(:,end-T_t:end)';
    series(2,:,:) = y_RE(:,end-T_t:end)';
    figname = [this_code, '_', 'sim' shocknames{s}, '_', gain_name, '_', PLM_name , '_', relevant_params, '_', date_today];
    subplot_names = titles_obs;
    legendnames = {'Learning', 'RE'};
    figtitle = [gain_title, '; shock imposed at t=', num2str(dt_vals(t))];
    create_subplot(series,subplot_names,figname,print_figs, figtitle, legendnames)
end

if plot_gains==1
    % 3) Average inverse gains
    yseries=mean(1./k,2)';
    xseries=1:T;
    seriesnames = '1/k';
    figname = [this_code, '_', 'loss','_', gain_name, '_', PLM_name ,  '_', relevant_params,'_', date_today];
    figtitle = ['Inverse gains ; ' , gain_title];
    create_plot(xseries,yseries,seriesnames,figname,print_figs,figtitle)
end