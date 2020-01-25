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

plot_simulated_sequence = 0;
skip_old_stuff = 1;

%% Parameters
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
p11 = param.p11;
p22 = param.p22;
ne = 3;

% Params for the general learning code
constant_only = 1; % learning constant only
mean_only_PLM = -1;
slope_and_constant = 2;

dgain = 1;  % 1 = decreasing gain, 2 = endogenous gain, 3 = constant gain
again_critCEMP  = 21;
again_critCUSUM = 22;
cgain = 3;

critCEMP=1;
critCUSUM=2;

%% Model selection and informational assumption

% Model selection
%%%%%%%%%%%%%%%%%%%
PLM = constant_only;
gain = again_critCEMP;
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

% Give names
if  PLM == constant_only
    PLM_name = 'constant_only'
elseif PLM == mean_only_PLM
    PLM_name = 'mean_only_PLM';
elseif PLM == slope_and_constant
    PLM_name = 'slope_and_constant'
end

if  gain == dgain
    gain_name = 'dgain'
    gain_title = 'Decreasing gain';
elseif gain == cgain
    gain_name = 'cgain'
    gain_title = 'Constant gain';
elseif gain == again_critCEMP
    gain_name = 'again_critCEMP';
    gain_title = 'Endogenous gain, CEMP''s criterion';
elseif gain == again_critCUSUM
    gain_name = 'again_critCUSUM';
    gain_title = 'Endogenous gain, CUSUM criterion';
end

%% Simulate models

% gen all the N sequences of shocks at once.
eN = randn(ne,T,N);

% Preallocate
nd = size(dt_vals,2);
d = 1; % the innovation, delta
GIR_Y_EE = zeros(ny,h,N,nd);
GIR_Y_LH = zeros(ny,h,N,nd);

% warning off
for s=2  %2->zoom in on monetary policy shock
    x0 = zeros(1,nx);
    x0(s) = d;
    
    for n=1:N
        % Sequence of innovationsﬂ
        e = squeeze(eN(:,:,n));
        
        % Unshocked
        dbstop if warning
        % RE
        [x_RE, y_RE] = sim_model(gx,hx,SIG,T,burnin,e);
        % Learning
        [x_LH, y_LH] = sim_learnLH(gx,hx,SIG,T,burnin,e, Aa, Ab, As, param, PLM, gain);
        
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
            [~, ys_LH,evening_fcst, morning_fcst, FA, FB, FEt_1, shock, diff] = sim_learnLH(gx,hx,SIG,T,burnin,e, Aa, Ab, As, param, PLM, gain, dt, x0);
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

if plot_simulated_sequence==1
clear series
T_t = 100;
% 2) SIMULATED HISTORY: OBSERVABLES LH against RE
series(1,:,:) = y_LH(:,end-T_t:end)';
series(2,:,:) = y_RE(:,end-T_t:end)';
figname = [this_code, '_', 'RIR_LH_' shocknames{s}, '_', gain_name, '_', PLM_name , '_', date_today];
subplot_names = titles_obs;
legendnames = {'Learning', 'RE'};
figtitle = [gain_title, '; shock imposed at t=', num2str(dt_vals(t))];
create_subplot(series,subplot_names,figname,print_figs, figtitle, legendnames)
end