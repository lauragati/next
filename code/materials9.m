% materials 8
% pretty much a copy of materials7 - i just want to be able to compare
% figures.
% Goals:
% 1. understand how to get i to go up on impact in RE
% 2. understand learning responses using expectation responses
% 1 Nov 2019
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
print_figs    = 0;
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
N = 10; % size of cross-section
dt_vals = [5,25]; % time of imposing innovation
d = 1; % the innovation, delta
Y_d = zeros(ny,T,N); % the unshocked observables  (ny x length of sim x size of cross-section)
Y_c = zeros(ny,T,N); 
F_d = zeros(4,T,N); % unshocked forecasts in this order: morning, evening, Fa, Fb
F_c = zeros(4,T,N);
YS_d = zeros(ny,T,N,2); % the shocked observables. The 4th dim reflects the different timing of imposing the shock.
YS_c = zeros(ny,T,N,2); 

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
        Y_d(:,:,n) =y_d;
        Y_c(:,:,n) =y_c;
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
            YS_d(:,:,n,t) =ys_d;
            YS_c(:,:,n,t) =ys_c;
        end
    end
end
disp(['(psi_x, psi_pi, rho, sig)=   ', num2str([psi_x, psi_pi, rho, sig])])
toc

return
%% GIRs for interest rate smoothing
d = 1; % the innovation, delta
shocknames = {'natrate', 'monpol','costpush'};
dt = 5; % when to impose shock
dt_vals = [5,25];
% cycle thru the shocks of the model
for s=2 %1:ne  %2->zoom in on monetary policy shock
    x0 = zeros(1,nx);
    x0(s) = d;
    h = 10; % horizon of IRF
    [IR, iry, irx]=ir(gx,hx,x0,h);
    iry = iry';
    RE_fcsts = gx*hx*irx';
    
    % For learning, the approach I take is resimulate everything and for each
    % period t, expose the econ to this same shock. Then take an average.
    GIRd = zeros(ny,h,T-h); % decreasing gain
    GIRc = zeros(ny,h,T-h); % constant gain
    GIR_fcst_d_m = zeros(h,T-h);
    GIR_fcst_c_m = zeros(h,T-h);
    GIR_fcst_d_e = zeros(h,T-h);
    GIR_fcst_c_e = zeros(h,T-h);
    for i=1:2
        t = dt_vals(i);
        % 1. create alternative simulations, adding the impulse always at either 5 or 25
        [~, ys_d, e_fcsts_d, m_fcsts_d, FA_ds, FB_ds] = sim_learn(gx,hx,SIG,T,burnin,e, Aa, Ab, As, param, setp, H, anal, constant_only, dgain, critCEMP,free, t, x0);
        [~, ys_c, e_fcsts_c, m_fcsts_c, FA_cs, FB_cs] = sim_learn(gx,hx,SIG,T,burnin,e, Aa, Ab, As, param, setp, H, anal, constant_only, cgain, critCEMP,free, t, x0);
        
        % 2. take differences between this and the standard simulation
        GIRd(:,:,t) = ys_d(:,t:t+h-1) - y_d(:,t:t+h-1);
        GIRc(:,:,t) = ys_c(:,t:t+h-1) - y_c(:,t:t+h-1);
        % do the same for expectations (just morning fcsts)
        GIR_fcst_d_m(:,t) = m_fcsts_d(t:t+h-1) - m_fcst_d(t:t+h-1);
        GIR_fcst_c_m(:,t) = m_fcsts_c(t:t+h-1) - m_fcst_c(t:t+h-1);
        GIR_fcst_d_e(:,t) = e_fcsts_d(t:t+h-1) - e_fcst_d(t:t+h-1); % % evening forecasts for decreasing learning
        GIR_fcst_c_e(:,t) = e_fcsts_c(t:t+h-1) - e_fcst_c(t:t+h-1); % evening forecasts for constant learning
        
        %LH expectations
        GIR_FA_d(:,:,t) = FA_ds(:,t:t+h-1) - FA_d(:,t:t+h-1);
        GIR_FB_d(:,:,t) = FB_ds(:,t:t+h-1) - FB_d(:,t:t+h-1);
        GIR_FA_c(:,:,t) = FA_cs(:,t:t+h-1) - FA_c(:,t:t+h-1);
        GIR_FB_c(:,:,t) = FB_cs(:,t:t+h-1) - FB_c(:,t:t+h-1);
    end
    
    
    % option 1: take simple averages
    RIRd = mean(GIRd,3);
    RIRc = mean(GIRc,3);
    
    RIRfa_d = mean(GIR_FA_d(1,:,:),3); % just take it for pi
    RIRfb_d = mean(GIR_FB_d(1,:,:),3);
    RIRfa_c = mean(GIR_FA_c(1,:,:),3);
    RIRfb_c = mean(GIR_FB_c(1,:,:),3);
    
    
    %     % option 2: sort and take percentile bands
    %     [lbd, medd, ubd] = confi_bands(GIRd,0.1);
    %     [lbc, medc, ubc] = confi_bands(GIRc,0.1);
    %
    %     % CI for forecasts
    %     [lb_mf_d, med_mf_d, ub_mf_d] = confi_bands(GIR_fcst_d,0.1);
    %     [lb_mf_c, med_mf_c, ub_mf_c] = confi_bands(GIR_fcst_c,0.1);
    %     [lb_ef_d, med_ef_d, ub_ef_d] = confi_bands(GIR_fcst_d_e,0.1);
    %     [lb_ef_c, med_ef_c, ub_ef_c] = confi_bands(GIR_fcst_c_e,0.1);
    
    return
    %% Plot IRFs
    titles = {'Inflation','Output gap','Int. rate', 'E^m_t(\pi_{t+1})', 'E^e_t(\pi_{t+1})'};
    
    % IRF observables constant gain learning only
    figure
    set(gcf,'color','w'); % sets white background color
    set(gcf, 'Position', get(0, 'Screensize')); % sets the figure fullscreen
    for i=1:ny
        subplot(1,ny,i)
        plot(zeros(1,h),'k--','linewidth', 2); hold on
        % Plot medians
        med = plot(medc(i,:),'color', dark_green,'linewidth', 2);
        % Plot CIs
        fill_Xcoord = [1:h, fliplr(1:h)];
        fillYcoord = [lbc(i,:), fliplr(ubc(i,:))];
        f = fill(fill_Xcoord, fillYcoord, light_green,'LineStyle','none');
        set(f,'facealpha',.5)
        title(titles(i))
        ax = gca; % current axes
        ax.FontSize = fs_prop;
        grid on
        grid minor
    end
    if print_figs ==1
        figname = [this_code, '_', 'IRFs_cgain_' shocknames{s},'_rho',rho_val, '_psi_pi_', psi_pi_val, '_sig_', sig_val, '_dt_', num2str(dt)]
        cd(figpath)
        export_fig(figname)
        cd(current_dir)
        close
    end
    
    % IRF observables decreasing gain learning only
    figure
    set(gcf,'color','w'); % sets white background color
    set(gcf, 'Position', get(0, 'Screensize')); % sets the figure fullscreen
    for i=1:ny
        subplot(1,ny,i)
        plot(zeros(1,h),'k--','linewidth', 2); hold on
        % Plot medians
        med = plot(medd(i,:),'b','linewidth', 2);
        % Plot CIs
        fill_Xcoord = [1:h, fliplr(1:h)];
        fillYcoord = [lbc(i,:), fliplr(ubc(i,:))];
        f = fill(fill_Xcoord, fillYcoord, light_green,'LineStyle','none');
        set(f,'facealpha',.5)
        title(titles(i))
        ax = gca; % current axes
        ax.FontSize = fs_prop;
        grid on
        grid minor
    end
    if print_figs ==1
        figname = [this_code, '_', 'IRFs_dgain_' shocknames{s},'_rho',rho_val, '_psi_pi_', psi_pi_val, '_sig_', sig_val, '_dt_', num2str(dt)]
        cd(figpath)
        export_fig(figname)
        cd(current_dir)
        close
    end
    
    % Plot evening and morning forecasts for constant gain only, IRF
    figure
    set(gcf,'color','w'); % sets white background color
    set(gcf, 'Position', get(0, 'Screensize')); % sets the figure fullscreen
    plot(zeros(1,h),'k--','linewidth', 2); hold on
    mf_c = plot(med_mf_c,'linewidth', 2);
    ef_c = plot(med_ef_c,'linewidth', 2);
    ylim([-0.15, 0])
    ax = gca; % current axes
    ax.FontSize = fs_prop;
    grid on
    grid minor
    legend([mf_c, ef_c],'morning fcst',  'evening fcst','location', 'southoutside')
    legend('boxoff')
    if print_figs ==1
        figname = [this_code, '_', 'IRFs_cgain_fcsts_me_' shocknames{s},'_rho',rho_val, '_psi_pi_', psi_pi_val, '_sig_', sig_val, '_dt_', num2str(dt)]
        cd(figpath)
        export_fig(figname)
        cd(current_dir)
        close
    end
    
    % Plot evening and morning forecasts for decreasing gain only, IRF
    figure
    set(gcf,'color','w'); % sets white background color
    set(gcf, 'Position', get(0, 'Screensize')); % sets the figure fullscreen
    plot(zeros(1,h),'k--','linewidth', 2); hold on
    mf_c = plot(med_mf_d,'linewidth', 2);
    ef_c = plot(med_ef_d,'linewidth', 2);
    ylim([-0.15, 0])
    ax = gca; % current axes
    ax.FontSize = fs_prop;
    grid on
    grid minor
    legend([mf_c, ef_c],'morning fcst',  'evening fcst','location', 'southoutside')
    legend('boxoff')
    if print_figs ==1
        figname = [this_code, '_', 'IRFs_dgain_fcsts_me_' shocknames{s},'_rho',rho_val, '_psi_pi_', psi_pi_val, '_sig_', sig_val, '_dt_', num2str(dt)]
        cd(figpath)
        export_fig(figname)
        cd(current_dir)
        close
    end
    
    % Plot LH expectations for constant and decreasing gain learning (just for pi)
    figure
    set(gcf,'color','w'); % sets white background color
    set(gcf, 'Position', get(0, 'Screensize')); % sets the figure fullscreen
    subplot(1,2,1)
    plot(zeros(1,h),'k--','linewidth', 2); hold on
    fa_d = plot(RIRfa_d,'linewidth', 2);
    fa_c = plot(RIRfa_c,'linewidth', 2);
    %     ylim([-0.15, 0])
    ax = gca; % current axes
    ax.FontSize = fs_prop;
    grid on
    grid minor
    legend([fa_d, fa_c],'decreasing',  'constant','location', 'southoutside')
    legend('boxoff')
    title('fa')
    subplot(1,2,2)
    plot(zeros(1,h),'k--','linewidth', 2); hold on
    fb_d = plot(RIRfb_d,'linewidth', 2);
    fb_c = plot(RIRfb_c,'linewidth', 2);
    %     ylim([-0.15, 0])
    ax = gca; % current axes
    ax.FontSize = fs_prop;
    grid on
    grid minor
    legend([fb_d, fb_c],'decreasing',  'constant','location', 'southoutside')
    legend('boxoff')
    title('fb')
    if print_figs ==1
        figname = [this_code, '_', 'IRFs_LH_' shocknames{s},'_rho',rho_val, '_psi_pi_', psi_pi_val, '_sig_', sig_val, '_dt_', num2str(dt)]
        cd(figpath)
        export_fig(figname)
        cd(current_dir)
        close
    end
    
end


