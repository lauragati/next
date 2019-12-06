% command_do_IRFs_for_cgain_dgain.m
% summarizes in a generic way the generation of IRFs for cgain and dgain LH
% learning.
% copies 1:1 materials10.m, but only does IRFs
% 21 Nov 2019


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

%% Parameters
tic
T = 400 % 400
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

% introduce adaptive names depending on the value of rho
rho_val_raw = num2str(rho);
rho_val = replace(rho_val_raw,'.','_');
psi_pi_val_raw = num2str(psi_pi);
psi_pi_val = replace(psi_pi_val_raw,'.','_');
rho_i_val_raw = num2str(rho_i);
rho_i_val = replace(rho_i_val_raw,'.','_');
alph_val_raw = num2str(alph);
alph_val = replace(alph_val_raw,'.','_');
sig_val_raw = num2str(sig);
sig_val = replace(sig_val_raw,'.','_');
gbar_val_raw = num2str(gbar);
gbar_val = replace(gbar_val_raw,'.','_');

current_param_values = [rho, rho_i, alph, kapp, psi_pi, sig, gbar]
current_param_names = ['\rho', '\rho_i', '\alpha', '\kappa', '\psi_{\pi}', '\sigma', '\bar{g}']

%% Model with interest rate smoothing

[fyn, fxn, fypn, fxpn] = model_NK_intrate_smoothing(param);
[gx,hx]=gx_hx_alt(fyn,fxn,fypn,fxpn);
[ny, nx] = size(gx);
[Ap_RE, As_RE, Aa, Ab, As] = matrices_A_intrate_smoothing(param, hx);


%% Simulate models w/ and w/o shocks for a cross-section

% Params for the general learning code
constant_only = 1; % learning constant only
mean_only_PLM = -1;
slope_and_constant = 2;
% lets alternate between these
PLM = constant_only;
if  PLM == constant_only
    PLM_name = 'constant_only';
elseif PLM == mean_only_PLM
    PLM_name = 'mean_only_PLM';
elseif PLM == slope_and_constant
    PLM_name = 'slope_and_constant';
end

dgain = 1;  % 1 = decreasing gain, 2 = endogenous gain, 3 = constant gain
again = 2;
cgain = 3;
critCEMP=1;
critCUSUM=2;

% Params specific to this cross-section exercise
N = 500; %500 size of cross-section
dt_vals = 25; % time of imposing innovation
nd = size(dt_vals,2);
d = 1; % the innovation, delta
h = 10; % h-period IRFs
nf = 4; % number of fcsts (morning, evening, fa and fb)
Y_d = zeros(ny,T,N); % the unshocked observables  (ny x length of sim x size of cross-section)
Y_c = zeros(ny,T,N);
F_d = zeros(nf,T,N); % unshocked forecasts in this order: morning, evening, Fa, Fb
F_c = zeros(nf,T,N);
FE_d = zeros(nf-2,T-1,N); % unshocked forecast errors, for pi only, 1st dim reflects morning and evening forecasts.
FE_c = zeros(nf-2,T-1,N);
FEL_d = zeros(T,N); % unshocked yesterday's evening FEs, for pi only.
FEL_c = zeros(T,N);
YS_d = zeros(ny,T,N,nd); % the shocked observables. The 4th dim reflects the different timing of imposing the shock.
YS_c = zeros(ny,T,N,nd);
FS_d = zeros(nf,T,N,nd); % the shocked forecasts in the same order; 4th dim is the two different impact times of shock.
FS_c = zeros(nf,T,N,nd);
FEs_d = zeros(nf-2,T-1,N,nd); % shocked forecast errors
FEs_c = zeros(nf-2,T-1,N,nd);
FELs_d = zeros(T,N,nd); % shocked yesterday's evening FEs, for pi only.
FELs_c = zeros(T,N,nd);
GIR_y_d = zeros(ny,h,N,nd);
GIR_y_c = zeros(ny,h,N,nd);
GIR_F_d = zeros(nf,h,N,nd);
GIR_F_c = zeros(nf,h,N,nd);
GIR_FE_d = zeros(nf-2,h,N,nd);
GIR_FE_c = zeros(nf-2,h,N,nd);
GIR_FEL_d = zeros(h,N,nd);
GIR_FEL_c = zeros(h,N,nd);

eN = randn(ne,T,N); % gen all the N sequences of shocks at once.
%     param = struct2array(setp);
%     fn = fieldnames(setp);
for s=2  %2->zoom in on monetary policy shock
    x0 = zeros(1,nx);
    x0(s) = d;
    h = 10; % horizon of IRF
    [IR, iry, irx]=ir(gx,hx,x0,h);
    iry = iry';
    RE_fcsts = gx*hx*irx';
    for n=1:N
        % Sequence of innovations
        e = [squeeze(eN(:,:,n)); zeros(1,T)]; % adding zero shocks to interest rate lag
        
        % Unshocked - let y denote unshocked
        [~, y_d, e_fcst_d, m_fcst_d, FA_d, FB_d, FEt_1_d] = sim_learn(gx,hx,SIG,T,burnin,e, Aa, Ab, As, param, PLM, dgain, critCEMP);
        [~, y_c, e_fcst_c, m_fcst_c, FA_c, FB_c, FEt_1_c] = sim_learn(gx,hx,SIG,T,burnin,e, Aa, Ab, As, param, PLM, cgain, critCEMP);
        
        % Gather unshocked observables - let big Y denote the the whole
        % cross-section, big F unshocked forecasts of inflation
        Y_d(:,:,n) = y_d;
        Y_c(:,:,n) = y_c;
        F_d(1,:,n) = m_fcst_d;
        F_d(2,:,n) = e_fcst_d;
        F_d(3,:,n) = FA_d(1,:);
        F_d(4,:,n) = FB_d(1,:);
        F_c(1,:,n) = m_fcst_c;
        F_c(2,:,n) = e_fcst_c;
        F_c(3,:,n) = FA_c(1,:);
        F_c(4,:,n) = FB_c(1,:);
        % construct unshocked forecast errors, for pi only
        FE_d(1,:,n) = Y_d(1,2:end,n) - F_d(1,1:end-1,n); % morning FE
        FE_d(2,:,n) = Y_d(1,2:end,n) - F_d(2,1:end-1,n); % evening FE
        FE_c(1,:,n) = Y_c(1,2:end,n) - F_c(1,1:end-1,n); % morning FE
        FE_c(2,:,n) = Y_c(1,2:end,n) - F_c(2,1:end-1,n); % evening FE
        % note: FEs go from t+1 ... T b/c the first one isn't defined since
        % they are realized at t+1.
        % Construct unshocked yesterday's evening FEs, for pi only
        FEL_d(:,n) = FEt_1_d;
        FEL_c(:,n) = FEt_1_c;
        
        for t=1:nd
            dt = dt_vals(t);
            % Shocked = let ys denote shocked
            [~, ys_d, e_fcsts_d, m_fcsts_d, FAs_d, FBs_d, FEst_1_d] = sim_learn(gx,hx,SIG,T,burnin,e, Aa, Ab, As, param, PLM, dgain, critCEMP, dt, x0);
            [~, ys_c, e_fcsts_c, m_fcsts_c, FAs_c, FBs_c, FEst_1_c] = sim_learn(gx,hx,SIG,T,burnin,e, Aa, Ab, As, param, PLM, cgain, critCEMP, dt, x0);
            
            % Gather shocked observables
            YS_d(:,:,n,t) = ys_d;
            YS_c(:,:,n,t) = ys_c;
            FS_d(1,:,n,t) = m_fcsts_d;
            FS_d(2,:,n,t) = e_fcsts_d;
            FS_d(3,:,n,t) = FAs_d(1,:);
            FS_d(4,:,n,t) = FBs_d(1,:);
            FS_c(1,:,n,t) = m_fcsts_c;
            FS_c(2,:,n,t) = e_fcsts_c;
            FS_c(3,:,n,t) = FAs_c(1,:);
            FS_c(4,:,n,t) = FBs_c(1,:);
            % construct shocked forecast errors, for pi only
            FEs_d(1,:,n,t) = YS_d(1,2:end,n,t) - FS_d(1,1:end-1,n,t);
            FEs_d(2,:,n,t) = YS_d(1,2:end,n,t) - FS_d(2,1:end-1,n,t);
            FEs_c(1,:,n,t) = YS_c(1,2:end,n,t) - FS_c(1,1:end-1,n,t);
            FEs_c(2,:,n,t) = YS_c(1,2:end,n,t) - FS_c(2,1:end-1,n,t);
            % Construct shocked yesterday's evening FEs, for pi only
            FELs_d(:,n,t) = FEst_1_d;
            FELs_c(:,n,t) = FEst_1_c;
            
            % Construct GIRs
            GIR_y_d(:,:,n,t) = YS_d(:,dt:dt+h-1,n,t) - Y_d(:,dt:dt+h-1,n);
            GIR_y_c(:,:,n,t) = YS_c(:,dt:dt+h-1,n,t) - Y_c(:,dt:dt+h-1,n);
            GIR_F_d(:,:,n,t) = FS_d(:,dt:dt+h-1,n,t) - F_d(:,dt:dt+h-1,n);
            GIR_F_c(:,:,n,t) = FS_c(:,dt:dt+h-1,n,t) - F_c(:,dt:dt+h-1,n);
            % think hard: FE_t realized at t+1
            GIR_FE_d(:,:,n,t) = FEs_d(:,dt:dt+h-1,n,t) - FE_d(:,dt:dt+h-1,n); % IRF of morning fcst starting period dt
            GIR_FE_c(:,:,n,t) = FEs_c(:,dt:dt+h-1,n,t) - FE_c(:,dt:dt+h-1,n); % IRF of evening fcst starting period dt
            % yesterday's FE_{t-1}, realized at t and used to update at t
            GIR_FEL_d(:,n,t) = FELs_d(dt:dt+h-1,n,t) - FEL_d(dt:dt+h-1,n);
            GIR_FEL_c(:,n,t) = FELs_c(dt:dt+h-1,n,t) - FEL_c(dt:dt+h-1,n);
        end
    end

% Construct RIRs by simple method: means (Option 1)
RIR_y_d = squeeze(mean(GIR_y_d,3));
RIR_y_c = squeeze(mean(GIR_y_c,3));
RIR_F_d = squeeze(nanmean(GIR_F_d,3));
RIR_F_c = squeeze(nanmean(GIR_F_c,3));
RIR_FE_d = squeeze(nanmean(GIR_FE_d,3));
RIR_FE_c = squeeze(nanmean(GIR_FE_c,3));
RIR_FEL_d = squeeze(nanmean(GIR_FEL_d,2));
RIR_FEL_c = squeeze(nanmean(GIR_FEL_c,2));

% % suggested value for gbar based on where decreasing gain is after 50 periods
% gbar_sugg = 1/dgain_at50;
% % Note that gbar_sugg - 50^(-1) (the std dgain algorithm) suggest very
% % close values

disp(['(psi_x, psi_pi, rho, sig)=   ', num2str([psi_x, psi_pi, rho, sig])])
toc

if stop_before_plots==1
    return
end
%% Plot 'em with new plot commands
shocknames = {'natrate', 'monpol','costpush'};
titles_obs = {'Inflation','Output gap','Int. rate'};
titles_fcsts = {'E^m_t(\pi_{t+1})', 'E^e_t(\pi_{t+1})'};
titles_LH = {'fa', 'fb'};
titles_FEs = {'FE^m_t(\pi_{t+1})', 'FE^e_t(\pi_{t+1})'};


    for t=1:nd % for the two diff times of imposing the shock
        dt = dt_vals(t);
        
        if skip_old_plots==0
            clear series
            % 1) OBSERVABLES DECREASING GAIN
            series(1,:,:) = RIR_y_d(:,:,t)';
            series(2,:,:) = iry';
            figname = [this_code, '_', 'RIR_y_d_' shocknames{s}, PLM_name,'_rho',rho_val, '_psi_pi_', psi_pi_val, '_sig_', sig_val, '_dt_', num2str(dt), '_gbar_', gbar_val];
            subplot_names = titles_obs;
            legendnames = {'Learning', 'RE'};
            figtitle = ['Decreasing gain, shock imposed at t=', num2str(dt_vals(t))];
            create_subplot(series,subplot_names,figname,print_figs, figtitle, legendnames)
            
            % 2) OBSERVABLES CONSTANT GAIN
            %             series = RIR_y_c(:,:,t);
            series(1,:,:) = RIR_y_c(:,:,t)';
            series(2,:,:) = iry';
            figname = [this_code, '_', 'RIR_y_c_' shocknames{s}, PLM_name,'_rho',rho_val, '_psi_pi_', psi_pi_val, '_sig_', sig_val, '_dt_', num2str(dt), '_gbar_', gbar_val];
            subplot_names = titles_obs;
            legendnames = {'Learning', 'RE'};
            figtitle = ['Constant gain, shock imposed at t=', num2str(dt_vals(t))];
            create_subplot(series,subplot_names,figname,print_figs, figtitle, legendnames)
            
            % 3) LH FORECASTS DECREASING GAIN
            series = RIR_F_d(3:end,:,t);
            figname = [this_code, '_', 'RIR_fafb_d_' shocknames{s}, PLM_name,'_rho',rho_val, '_psi_pi_', psi_pi_val, '_sig_', sig_val, '_dt_', num2str(dt), '_gbar_', gbar_val];
            subplot_names = titles_LH;
            figtitle = ['Decreasing gain, shock imposed at t=', num2str(dt_vals(t))];
            create_subplot(series,subplot_names,figname,print_figs, figtitle)
            
            % 4) LH FORECASTS CONSTANT GAIN
            series = RIR_F_c(3:end,:,t);
            figname = [this_code, '_', 'RIR_fafb_c_' shocknames{s}, PLM_name,'_rho',rho_val, '_psi_pi_', psi_pi_val, '_sig_', sig_val, '_dt_', num2str(dt), '_gbar_', gbar_val];
            subplot_names = titles_LH;
            figtitle = ['Constant gain, shock imposed at t=', num2str(dt_vals(t))];
            create_subplot(series,subplot_names,figname,print_figs, figtitle)
        end % skip old plots
        
        % 5) FORECASTS
        subplot1 = [RIR_F_d(1,:,t); RIR_F_d(2,:,t)];
        subplot2 = [RIR_F_c(1,:,t); RIR_F_c(2,:,t)];
        series = cat(3,subplot1, subplot2);
        figname = [this_code, '_', 'RIR_F_both_' shocknames{s}, PLM_name,'_rho',rho_val, '_psi_pi_', psi_pi_val, '_sig_', sig_val, '_dt_', num2str(dt), '_gbar_', gbar_val];
        subplot_names = {'Decreasing gain', 'Constant gain'};
        legend_entries = {'Morning', 'Evening'};
        figtitle = ['1-period ahead fcsts, shock imposed at t=', num2str(dt_vals(t))];
        create_subplot(series,subplot_names,figname,print_figs,figtitle,legend_entries)
        
        
        % 6) FORECAST ERRORS
        subplot1 = [RIR_FE_d(1,:,t); RIR_FE_d(2,:,t);RIR_FEL_d(:,t)'];
        subplot2 = [RIR_FE_c(1,:,t); RIR_FE_c(2,:,t);RIR_FEL_c(:,t)'];
        series = cat(3,subplot1, subplot2);
        figname = [this_code, '_', 'RIR_FE_both_' shocknames{s}, PLM_name,'_rho',rho_val, '_psi_pi_', psi_pi_val, '_sig_', sig_val, '_dt_', num2str(dt), '_gbar_', gbar_val];
        subplot_names = {'Decreasing gain', 'Constant gain'};
        legend_entries = {'Morning', 'Evening', 'Yesterday evening'};
        figtitle = ['1-period ahead FEs, shock imposed at t=', num2str(dt_vals(t))];
        create_subplot(series,subplot_names,figname,print_figs,figtitle,legend_entries)
        
    end
end



