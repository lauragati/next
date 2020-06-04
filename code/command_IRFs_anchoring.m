% command_IRFs_anchoring
% Do IRFs for anchoring true-baseline learning models
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

plot_IRFs=1;
plot_simulated_sequence = 1;
plot_gains=1;
plot_gain_IRF = 0;
plot_IRFs_anch = 0;
skip_old_stuff = 1;

%% Parameters
tic
[param, set, param_names, param_values_str, param_titles] = parameters_next;


psi_x = param.psi_x;
psi_pi = param.psi_pi;
thetbar =param.thetbar;
thettilde = param.thettilde;

sig_r = param.sig_r;
sig_i = param.sig_i;
sig_u = param.sig_u;
ne = 3;

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


%% Model selection and informational assumption

% Model selection
%%%%%%%%%%%%%%%%%%%
PLM = constant_only;
gain = again_critCUSUM;
%%%%%%%%%%%%%%%%%%%

T = 400 % 400
% Size of cross-section
N = 100 %100, 500
burnin = 0; % 100
dt_vals = 25; %25 time of imposing innovation 345
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
eN = randn(ne,T+burnin,N);

% Preallocate
nd = size(dt_vals,2);
d = 1; % the innovation, delta
GIR_Y_LH = zeros(ny,h,N,nd);
GIR_k = zeros(h,N,nd);
k = zeros(T,N);
ks= zeros(T,N);
anch = zeros(1,N); % indicator for whether that sequence was anchored when the shock hit

% warning off
for s=2  %2->zoom in on monetary policy shock
    x0 = zeros(1,nx);
    x0(s) = d;
    
    for n=1:N
        % Sequence of innovations
        e = squeeze(eN(:,:,n));
        
        % Unshocked
        %         dbstop if warning
        % RE
        [x_RE, y_RE] = sim_model(gx,hx,SIG,T,burnin,e);
        % Learning
%         dbstop in sim_learnLH at 116 if t>=4
        [x_LH, y_LH, ~, ~, ~, ~, ~, ~, diff,~, k(:,n)] = sim_learnLH(gx,hx,SIG,T+burnin,burnin,e, Aa, Ab, As, param, PLM, gain);
        
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
            [~, ys_LH, ~, ~, ~, ~, ~, ~, ~,~, ks(:,n), anch(n)] = sim_learnLH(gx,hx,SIG,T+burnin,burnin,e, Aa, Ab, As, param, PLM, gain, dt, x0);
            % Construct GIRs
            GIR_Y_LH(:,:,n,t) = ys_LH(:,dt:dt+h-1) - y_LH(:,dt:dt+h-1);
            GIR_k(:,n,t) = ks(dt:dt+h-1,n) - k(dt:dt+h-1,n);
        end
        
    end
end
% warning on
% Construct RIRs by simple method: means (Option 1)
RIR_Y_LH = squeeze(mean(GIR_Y_LH,3));
RIR_anch = squeeze(mean(GIR_Y_LH(:,:,find(anch)),3));
RIR_unanch = squeeze(mean(GIR_Y_LH(:,:,find(anch-1)),3));
RIR_k = squeeze(mean(GIR_k,2));
RIR_kinv = RIR_k;
% only invert for nonzero elements
RIR_kinv(abs(RIR_k) > 0) = 1./RIR_k(abs(RIR_k) > 0);

disp(['(psi_x, psi_pi, thetbar, thettilde)=   ', num2str([psi_x, psi_pi, thetbar, thettilde])])
toc

if stop_before_plots==1
    return
end


%% Plots

shocknames = {'natrate', 'monpol','costpush'};
titles_obs = {'Inflation','Output gap','Interest rate'};
titles_fcsts = {'E^m_t(\pi_{t+1})', 'E^e_t(\pi_{t+1})'};
titles_LH = {'fa', 'fb'};
titles_FEs = {'FE^m_t(\pi_{t+1})', 'FE^e_t(\pi_{t+1})'};

% create indexes for the positions of the parameters
fn = fieldnames(param);
make_index(fn)
interesting_param_names = param_names([psi_pi_idx, psi_x_idx, gbar_idx, thetbar_idx, thettilde_idx, kap_idx, lamx_idx, lami_idx]);
interesting_param_vals = param_values_str([psi_pi_idx, psi_x_idx, gbar_idx, thetbar_idx, thettilde_idx, kap_idx,  lamx_idx, lami_idx]);
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
        legendnames = {'Anchoring', 'RE'};
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
    legendnames = {'Anchoring', 'RE'};
    figtitle = [gain_title, '; shock imposed at t=', num2str(dt_vals(t))];
    create_subplot(series,subplot_names,figname,print_figs, figtitle, legendnames)
end

if plot_gains==1
    % 3) Average inverse gains
    yseries=mean(1./k,2)';
    xseries=1:T;
    seriesnames = 'k^{-1}';
    figname = [this_code, '_', 'loss','_', gain_name, '_', PLM_name ,  '_', relevant_params,'_', date_today];
    figtitle = ['Gains ; ' , gain_title];
    %     figtitle = '';
    create_plot(xseries,yseries,seriesnames,figname,print_figs,figtitle)
end

if plot_gain_IRF==1
    % 4) Inverse gain IRF
    clear yseries xseries
    yseries=RIR_kinv';
    xseries=1:h;
    seriesnames = '1/k';
    figname = [this_code, '_', 'loss','_', gain_name, '_', PLM_name ,  '_', relevant_params,'_', date_today];
    figtitle = ['IRF gains ; ' , gain_title];
    create_plot(xseries,yseries,seriesnames,figname,print_figs,figtitle)
end

if plot_IRFs_anch==1
    for t=1:nd % for the two diff times of imposing the shock
        dt = dt_vals(t);
        
        clear series
        % 5) IRF: OBSERVABLES LH against RE, anchored
        series(1,:,:) = RIR_anch(:,:,t)';
        series(2,:,:) = iry';
        figname = [this_code, '_', 'RIR_LH_anch_' shocknames{s}, '_', gain_name, '_', PLM_name , '_', date_today];
        subplot_names = titles_obs;
        legendnames = {'Learning', 'RE'};
        figtitle = [gain_title, '; when shock imposed at t=', num2str(dt_vals(t)), ', anchored'];
%         figtitle = '';
        create_subplot(series,subplot_names,figname,print_figs, figtitle, legendnames)
        
        clear series
        % 6) IRF: OBSERVABLES LH against RE, unanchored
        series(1,:,:) = RIR_unanch(:,:,t)';
        series(2,:,:) = iry';
        figname = [this_code, '_', 'RIR_LH_unanch_' shocknames{s}, '_', gain_name, '_', PLM_name , '_', date_today];
%         figname = [this_code, '_', 'RIR_LH_unanch_' shocknames{s}, '_', gain_name, '_', PLM_name , '_', param_names_vals{1}, '_', date_today];
        subplot_names = titles_obs;
        legendnames = {'Anchoring', 'RE'};
        figtitle = [gain_title, '; when shock imposed at t=', num2str(dt_vals(t)), ', unanchored'];
%         figtitle = '';
        create_subplot(series,subplot_names,figname,print_figs, figtitle, legendnames)
    end
end

figure
plot(diff)
title('Convergence')