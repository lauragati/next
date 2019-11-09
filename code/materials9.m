% materials 8
% Goals:
% 1. understand overshooting better - constant vs decreasing gain
% 8 Nov 2019
clearvars
close all
clc

% Add all the relevant paths
current_dir = pwd;
cd ../ % go up 1 levels
basepath = pwd;
cd .. % go up another level to BC_Research
BC_researchpath = pwd;
toolpath = [BC_researchpath '/matlab_toolbox'];
export_figpath = [toolpath '/Export_Fig'];
figpath = [basepath '/figures'];
tablepath = [basepath '/tables'];
datapath = [basepath '/data'];

cd(current_dir)

addpath(basepath)
addpath(toolpath)
addpath(export_figpath)
addpath(figpath)
addpath(datapath)

this_code = mfilename;

% Variable stuff ---
print_figs    = 1;
do_old_plots  = 0;
if print_figs ==1
    output_table  = 1;
else
    output_table =0;
end
skip_old_plots =1;

fs=20; % fontsize
lw=2; % linewidth
fs_pres = 80;
lw_pres = 6;
fs_prop = 40;
lw_prop = 4;
% Some color spectra
% grey color (divide by 255)
grey = [128,128,128]/255;
silver = [192,192,192]/255;
maroon = [138 0 0]/255;
light_coral = [240 128 128]/255;
light_salmon = [255,160,122]/255;
dark_green = [0 100 0]/255;
green = [0 128 0]/255;
light_green = [144 238 144]/255;
light_sky_blue = [135 206 250]/255;

teal = [0,128,128]/255;
purple = [128,0,128]/255;
saddle_brown = [139,69,19]/255;
light_brown = [181,101,29]/255;

%% Simulation
tic
T = 400
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


%% Model with interest rate smoothing

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

current_param_values = [rho, rho_i, alph, kapp, psi_pi, sig]
current_param_names = ['\rho', '\rho_i', '\alpha', '\kappa', '\psi_{\pi}', '\sigma']

[fyn, fxn, fypn, fxpn] = model_NK_intrate_smoothing(param);
[gx,hx]=gx_hx_alt(fyn,fxn,fypn,fxpn);
[ny, nx] = size(gx);
[Ap_RE, As_RE, Aa, Ab, As] = matrices_A_intrate_smoothing(param, setp, hx);


%% Simulate models w/ and w/o shocks for a cross-section

% Simulate models
% Use Ryan's code to simulate from the RE model
% [x_RE, y_RE] = sim_model(gx,hx,SIG,T,burnin,e);

% Params for the general learning code
anal = 1; % take analytical LR exp
H=0;
gain = 2; % 1 = decreasing gain, 2 = endogenous gain, 3 = constant gain
criterion = 2; % 1= CEMP's, 2 = CUSUM
% dt = 0; % when shock imposed. If zero or not specified, then no shock
constant_only = 1; % learning constant only
dgain = 1;
again = 2;
cgain = 3;
critCEMP=1;
critCUSUM=2;
free=1; % use versions of the code that are n-free (use hx instead of P)
not_free=0;

% Params specific to this cross-section exercise
N = 500; % size of cross-section
dt_vals = [5,25]; % time of imposing innovation
d = 1; % the innovation, delta
h = 10; % h-period IRFs
nf = 4; % number of fcsts (morning, evening, fa and fb)
Y_d = zeros(ny,T,N); % the unshocked observables  (ny x length of sim x size of cross-section)
Y_c = zeros(ny,T,N);
F_d = zeros(nf,T,N); % unshocked forecasts in this order: morning, evening, Fa, Fb
F_c = zeros(nf,T,N);
YS_d = zeros(ny,T,N,2); % the shocked observables. The 4th dim reflects the different timing of imposing the shock.
YS_c = zeros(ny,T,N,2);
FS_d = zeros(nf,T,N,2); % the shocked forecasts in the same order; 4th dim is the two different impact times of shock.
FS_c = zeros(nf,T,N,2);
GIR_y_d = zeros(ny,h,N,2);
GIR_y_c = zeros(ny,h,N,2);
GIR_F_d = zeros(nf,h,N,2);
GIR_F_c = zeros(nf,h,N,2);

for s=2 %1:ne  %2->zoom in on monetary policy shock
    x0 = zeros(1,nx);
    x0(s) = d;
    h = 10; % horizon of IRF
    [IR, iry, irx]=ir(gx,hx,x0,h);
    iry = iry';
    RE_fcsts = gx*hx*irx';
    for n=1:N
        % Sequence of innovations
        e = randn(ne,T);
        e = [e; zeros(1,T)]; % adding zero shocks to interest rate lag
        
        % Unshocked - let y denote unshocked
        [~, y_d, e_fcst_d, m_fcst_d, FA_d, FB_d] = sim_learn(gx,hx,SIG,T,burnin,e, Aa, Ab, As, param, setp, H, anal, constant_only, dgain, critCEMP, free);
        [~, y_c, e_fcst_c, m_fcst_c, FA_c, FB_c] = sim_learn(gx,hx,SIG,T,burnin,e, Aa, Ab, As, param, setp, H, anal, constant_only, cgain, critCEMP, free);
        
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
        
        for t=1:2
            dt = dt_vals(t);
            % Shocked = let ys denote shocked
            [~, ys_d, e_fcsts_d, m_fcsts_d, FAs_d, FBs_d] = sim_learn(gx,hx,SIG,T,burnin,e, Aa, Ab, As, param, setp, H, anal, constant_only, dgain, critCEMP,free, dt, x0);
            [~, ys_c, e_fcsts_c, m_fcsts_c, FAs_c, FBs_c] = sim_learn(gx,hx,SIG,T,burnin,e, Aa, Ab, As, param, setp, H, anal, constant_only, cgain, critCEMP,free, dt, x0);
            
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
            
            % Construct GIRs
            GIR_y_d(:,:,n,t) = YS_d(:,dt:dt+h-1,n,t) - Y_d(:,dt:dt+h-1,n);
            GIR_y_c(:,:,n,t) = YS_c(:,dt:dt+h-1,n,t) - Y_c(:,dt:dt+h-1,n);
            GIR_F_d(:,:,n,t) = FS_d(:,dt:dt+h-1,n,t) - F_d(:,dt:dt+h-1,n);
            GIR_F_c(:,:,n,t) = FS_c(:,dt:dt+h-1,n,t) - F_c(:,dt:dt+h-1,n);
        end
    end
end

% Construct RIRs by simple method: means (Option 1)
RIR_y_d = squeeze(mean(GIR_y_d,3));
RIR_y_c = squeeze(mean(GIR_y_c,3));
RIR_F_d = squeeze(mean(GIR_F_d,3));
RIR_F_c = squeeze(mean(GIR_F_c,3));

disp(['(psi_x, psi_pi, rho, sig)=   ', num2str([psi_x, psi_pi, rho, sig])])
toc


%% Plot 'em
shocknames = {'natrate', 'monpol','costpush'};
titles_obs = {'Inflation','Output gap','Int. rate', 'E^m_t(\pi_{t+1})'};
titles_fcsts = {'E^m_t(\pi_{t+1})', 'E^e_t(\pi_{t+1})', 'fa', 'fb'};

for s=2 % for all the shocks (right now only monpol)
    for t=1:2 % for the two diff times of imposing the shock
        dt = dt_vals(t);
        % 1) Observables for decreasing gain
        figure
        set(gcf,'color','w'); % sets white background color
        set(gcf, 'Position', get(0, 'Screensize')); % sets the figure fullscreen
        sgtitle(['Decreasing gain, shock imposed at t=', num2str(dt_vals(t))], 'FontSize',fs_prop )
        for i=1:ny
            subplot(1,ny,i)
            plot(zeros(1,h),'k--','linewidth', 2); hold on
            h1 = plot(RIR_y_d(i,:,t),'linewidth', 2);
            ax = gca; % current axes
            ax.FontSize = fs_prop;
            grid on
            grid minor
            title(titles_obs(i))
        end
        if print_figs ==1
            figname = [this_code, '_', 'RIR_y_d_' shocknames{s},'_rho',rho_val, '_psi_pi_', psi_pi_val, '_sig_', sig_val, '_dt_', num2str(dt)]
            cd(figpath)
            export_fig(figname)
            cd(current_dir)
            close
        end
        
        % 2.) Observables for constant gain
        figure
        set(gcf,'color','w'); % sets white background color
        set(gcf, 'Position', get(0, 'Screensize')); % sets the figure fullscreen
        sgtitle(['Constant gain, shock imposed at t=', num2str(dt_vals(t))], 'FontSize',fs_prop )
        for i=1:ny
            subplot(1,ny,i)
            plot(zeros(1,h),'k--','linewidth', 2); hold on
            h1 = plot(RIR_y_c(i,:,t),'linewidth', 2);
            ax = gca; % current axes
            ax.FontSize = fs_prop;
            grid on
            grid minor
            title(titles_obs(i))
        end
        if print_figs ==1
            figname = [this_code, '_', 'RIR_y_c_' shocknames{s},'_rho',rho_val, '_psi_pi_', psi_pi_val, '_sig_', sig_val, '_dt_', num2str(dt)]
            cd(figpath)
            export_fig(figname)
            cd(current_dir)
            close
        end
        
        % Forecasts for decreasing gain
        figure
        set(gcf,'color','w'); % sets white background color
        set(gcf, 'Position', get(0, 'Screensize')); % sets the figure fullscreen
        sgtitle(['Decreasing gain, shock imposed at t=', num2str(dt_vals(t))], 'FontSize',fs_prop )
        for i=1:nf
            subplot(1,nf,i)
            plot(zeros(1,h),'k--','linewidth', 2); hold on
            h1 = plot(RIR_F_d(i,:,t),'linewidth', 2);
            ax = gca; % current axes
            ax.FontSize = fs_prop;
            grid on
            grid minor
            title(titles_fcsts(i))
        end
        if print_figs ==1
            figname = [this_code, '_', 'RIR_F_d_' shocknames{s},'_rho',rho_val, '_psi_pi_', psi_pi_val, '_sig_', sig_val, '_dt_', num2str(dt)]
            cd(figpath)
            export_fig(figname)
            cd(current_dir)
            close
        end
        
        % Forecasts for decreasing gain
        figure
        set(gcf,'color','w'); % sets white background color
        set(gcf, 'Position', get(0, 'Screensize')); % sets the figure fullscreen
        sgtitle(['Constant gain, shock imposed at t=', num2str(dt_vals(t))], 'FontSize',fs_prop )
        for i=1:nf
            subplot(1,nf,i)
            plot(zeros(1,h),'k--','linewidth', 2); hold on
            h1 = plot(RIR_F_c(i,:,t),'linewidth', 2);
            ax = gca; % current axes
            ax.FontSize = fs_prop;
            grid on
            grid minor
            title(titles_fcsts(i))
        end
        if print_figs ==1
            figname = [this_code, '_', 'RIR_F_c_' shocknames{s},'_rho',rho_val, '_psi_pi_', psi_pi_val, '_sig_', sig_val, '_dt_', num2str(dt)]
            cd(figpath)
            export_fig(figname)
            cd(current_dir)
            close
        end
        
    end
end
