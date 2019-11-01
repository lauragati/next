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

% Sequence of innovations
rng(0)
e = randn(ne,T);
e = [e; zeros(1,T)]; % adding zero shocks to interest rate lag


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

current_param_values = [rho, rho_i, alph, kapp, psi_pi]
current_param_names = ['\rho', '\rho_i', '\alpha', '\kappa', '\psi_{\pi}']
% Solve RE model
% First old version w/o int rate smoothing
[fyn, fxn, fypn, fxpn] = model_NK(param);
[gx0,hx0]=gx_hx_alt(fyn,fxn,fypn,fxpn);

[fyn, fxn, fypn, fxpn] = model_NK_intrate_smoothing(param);
[gx,hx]=gx_hx_alt(fyn,fxn,fypn,fxpn);
[ny, nx] = size(gx);
[Ap_RE, As_RE, Aa, Ab, As] = matrices_A_intrate_smoothing(param, setp, hx); % perfect - Aa, Ab are the same as before, As has extra column of zeros, otherwise identical!
bigguy = -[fyn fxpn]^(-1)*[fypn fxn];
% Simulate models
% Use Ryan's code to simulate from the RE model
[x_RE, y_RE] = sim_model(gx,hx,SIG,T,burnin,e);


% Now simulate the learning models using the general code w/o imposing
% shocks
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
[x_d, y_d, e_fcst_d, m_fcst_d] = sim_learn(gx,hx,SIG,T,burnin,e, Aa, Ab, As, param, setp, H, anal, constant_only, dgain, critCEMP, free);
[x_a, y_a, e_fcst_a, m_fcst_a, ~, ~,pibar_a, k_a] = sim_learn(gx,hx,SIG,T,burnin,e, Aa, Ab, As, param, setp, H, anal, constant_only, again, critCEMP, free);
[x_a_cusum, y_a_cusum,e_fcst_ac, m_fcst_ac, ~, ~,pibar_a_cusum, k_cusum] = sim_learn(gx,hx,SIG,T,burnin,e, Aa, Ab, As, param, setp, H, anal, constant_only, again, critCUSUM, free);
[x_c, y_c, e_fcst_c, m_fcst_c] = sim_learn(gx,hx,SIG,T,burnin,e, Aa, Ab, As, param, setp, H, anal, constant_only, cgain, critCEMP, free);
ka_inv = 1./k_a;
k_cusum_inv = 1./k_cusum;


% Gather observables
Z(:,:,1) =y_RE;
Z(:,:,2) =y_d;
Z(:,:,3) =y_a;
Z(:,:,4) =y_c;
Z(:,:,5) =y_a_cusum;


%% GIRs for interest rate smoothing
d = 1; % the innovation, delta
shocknames = {'natrate', 'monpol','costpush'};
titles = {'Inflation','Output gap','Int. rate', 'E^m_t(\pi_{t+1})'};

% % Let's check shocks in RE
% for s=1:ne
%     x0 = zeros(1,nx);
%     x0(s) = d;
%     h = 10; % horizon of IRF
%     [IR, iry, irx]=ir(gx,hx,x0,h);
%     iry = iry';
%
%     figure
%     set(gcf,'color','w'); % sets white background color
%     set(gcf, 'Position', get(0, 'Screensize')); % sets the figure fullscreen
%     subplot(1,3,1)
%     plot(zeros(1,h),'k--','linewidth', lw_prop); hold on
%     plot(iry(1,:),'linewidth', lw_prop)
%     title(titles(1))
%     ax = gca; % current axes
%     ax.FontSize = fs_prop;
%     grid on
%     grid minor
%     subplot(1,3,2)
%     plot(zeros(1,h),'k--','linewidth', lw_prop); hold on
%     plot(iry(2,:),'linewidth', lw_prop)
%     title(titles(2))
%     ax = gca; % current axes
%     ax.FontSize = fs_prop;
%     grid on
%     grid minor
%     subplot(1,3,3)
%     plot(zeros(1,h),'k--','linewidth', lw_prop); hold on
%     plot(iry(3,:),'linewidth', lw_prop)
%     title(titles(3))
%     ax = gca; % current axes
%     ax.FontSize = fs_prop;
%     grid on
%     grid minor
%     if print_figs ==1
%         figname = [this_code, '_', 'IRFs_RE_' shocknames{s},'_rho_',rho_val, '_rho_i_', rho_i_val, '_alpha_', alph_val]
%         cd(figpath)
%         export_fig(figname)
%         cd(current_dir)
%         close
%     end
% end
%
% current_param_names
% current_param_values
%
% return
% cycle thru the shocks of the model
for s=2 %1:ne zoom in on monetary policy shock 
    x0 = zeros(1,nx);
    x0(s) = d;
    h = 10; % horizon of IRF
    [IR, iry, irx]=ir(gx,hx,x0,h);
    iry = iry';
    RE_fcsts = gx*hx*irx';
    
    % For learning, the approach I take is resimulate everything and for each
    % period t, expose the econ to this same shock. Then take an average.
    GIRd = zeros(ny,h,T-h); % decreasing gain
    GIRa = zeros(ny,h,T-h); % anchoring
    GIRc = zeros(ny,h,T-h); % constant gain
    pibars_a = zeros(T,T-h);
    pibars_a_cusum = zeros(T,T-h);
    ks_a = zeros(T,T-h);
    ks_a_cusum = zeros(T,T-h);
    GIR_fcst_d = zeros(h,T-h);
    GIR_fcst_a = zeros(h,T-h);
    GIR_fcst_c = zeros(h,T-h);
    GIR_k = zeros(h,T-h);
    GIR_k_cusum = zeros(h,T-h);
    GIR_k2 = zeros(h,T-h);
    for t=2:T-h
        % 1. create alternative simulations, adding the impulse always at a new time t
        % Now let's shock our general learning code
        [xs_d, ys_d, e_fcsts_d, m_fcsts_d] = sim_learn(gx,hx,SIG,T,burnin,e, Aa, Ab, As, param, setp, H, anal, constant_only, dgain, critCEMP,free, t, x0);
        [xs_a, ys_a, e_fcsts_a, m_fcsts_a, ~, ~,pibars_a(:,t), ks_a(:,t)] = sim_learn(gx,hx,SIG,T,burnin,e, Aa, Ab, As, param, setp, H, anal, constant_only, again, critCEMP,free, t, x0);
        [xs_c, ys_c, e_fcsts_c, m_fcsts_c] = sim_learn(gx,hx,SIG,T,burnin,e, Aa, Ab, As, param, setp, H, anal, constant_only, cgain, critCEMP,free, t, x0);
        [xs_a_cusum, ys_a_cusum,~,~, ~, ~,pibars_a_cusum(:,t), ks_a_cusum(:,t)] = sim_learn(gx,hx,SIG,T,burnin,e, Aa, Ab, As, param, setp, H, anal, constant_only, again, critCUSUM,free, t, x0);
        
        % 2. take differences between this and the standard simulation
        GIRd(:,:,t) = ys_d(:,t:t+h-1) - y_d(:,t:t+h-1);
        GIRa(:,:,t) = ys_a(:,t:t+h-1) - y_a(:,t:t+h-1);
        GIRc(:,:,t) = ys_c(:,t:t+h-1) - y_c(:,t:t+h-1);
        % do the same for expectations (just morning fcsts)
        GIR_fcst_d(:,t) = m_fcsts_d(t:t+h-1) - m_fcst_d(t:t+h-1);
        GIR_fcst_a(:,t) = m_fcsts_a(t:t+h-1) - m_fcst_a(t:t+h-1);
        GIR_fcst_c(:,t) = m_fcsts_c(t:t+h-1) - m_fcst_c(t:t+h-1);
        % and for the gains
        GIR_k(:,t) = 1./ks_a(t:t+h-1,t)' - 1./k_a(t:t+h-1);
        GIR_k_cusum(:,t) = 1./ks_a_cusum(t:t+h-1,t)' - 1./k_cusum(t:t+h-1);
        
        GIR_k2(:,t) = ks_a(t:t+h-1,t)' - k_a(t:t+h-1);
    end
    
    % option 1: take simple averages
    RIRd = mean(GIRd,3);
    RIRa = mean(GIRa,3);
    RIRc = mean(GIRc,3);
    
    RIRk = mean(GIR_k,2);
    RIRk_cusum = mean(GIR_k_cusum,2);
    
    % take averages of gains and drifts too
    pibars_a_mean = mean(pibars_a,2);
    pibars_a_cusum_mean = mean(pibars_a_cusum,2);
    ks_a_mean = mean(ks_a,2);
    ks_a_cusum_mean = mean(ks_a_cusum,2);
    inv_ks_a_mean = 1./ks_a_mean;
    inv_ks_a_cusum_mean = 1./ks_a_cusum_mean;
    
    % option 2: sort and take percentile bands
    [lbd, medd, ubd] = confi_bands(GIRd,0.1);
    [lba, meda, uba] = confi_bands(GIRa,0.1);
    [lbc, medc, ubc] = confi_bands(GIRc,0.1);
    
    % CI for k doesn't make a lot of sense because of the inverse which
    % makes everything pretty darn small 
    [lb_k, med_k, ub_k] = confi_bands(1./ks_a,0.1);
    [lb_k_cusum, med_k_cusum, ub_k_cusum] = confi_bands(1./ks_a_cusum,0.1);
    
    % ... or  you explode
    [lb_k2, med_k2, ub_k2] = confi_bands(GIR_k2,0.01);
    inv_lb_k2 = 1./lb_k2;
    inv_ub_k2 = 1./ub_k2;
    
    % CI for forecasts
    [lb_mf_d, med_mf_d, ub_mf_d] = confi_bands(GIR_fcst_d,0.1);
    [lb_mf_a, med_mf_a, ub_mf_a] = confi_bands(GIR_fcst_a,0.1);
    [lb_mf_c, med_mf_c, ub_mf_c] = confi_bands(GIR_fcst_c,0.1);
    
    %% Plot IRFs
    titles = {'Inflation','Output gap','Int. rate', 'E^m_t(\pi_{t+1})'};
    
    figure
    set(gcf,'color','w'); % sets white background color
    set(gcf, 'Position', get(0, 'Screensize')); % sets the figure fullscreen
    for i=1:ny+1
        if i>ny
            subplot(1,ny+1,i)
            plot(zeros(1,h),'k--','linewidth', 2); hold on
            fcst_RE = plot(RE_fcsts(1,:),'k','linewidth', 2);
            mf_d = plot(med_mf_d,'b','linewidth', 2);
            mf_a = plot(med_mf_a,'r','linewidth', 2);
            mf_c = plot(med_mf_c,'color', dark_green,'linewidth', 2);
            legend([fcst_RE, mf_d, mf_a, mf_c],'RE', 'Decreasing', 'Anchor', 'Constant ','location', 'southoutside')
            legend('boxoff')
        else
            subplot(1,ny+1,i)
            plot(zeros(1,h),'k--','linewidth', 2); hold on
            re = plot(iry(i,:),'k','linewidth', 2);
            % Plot means
            %         d_mean = plot(RIRd(i,:),'b','linewidth', 2);
            %         am_mean = plot(RIRa(i,:),'r','linewidth', 2);
            %         c_mean = plot(RIRc(i,:),'color', dark_green,'linewidth', 2);
            
            % Plot medians
            d_med = plot(medd(i,:),'b','linewidth', 2);
            a_med = plot(meda(i,:),'r','linewidth', 2);
            c_med = plot(medc(i,:),'color', dark_green,'linewidth', 2);
            % Plot CIs
            fill_Xcoord = [1:h, fliplr(1:h)];
            fillYcoord = [lbd(i,:), fliplr(ubd(i,:))];
            f = fill(fill_Xcoord, fillYcoord, light_sky_blue,'LineStyle','none');
            set(f,'facealpha',.5)
            fillYcoord = [lba(i,:), fliplr(uba(i,:))];
            f = fill(fill_Xcoord, fillYcoord, light_salmon,'LineStyle','none');
            set(f,'facealpha',.5)
            fillYcoord = [lbc(i,:), fliplr(ubc(i,:))];
            f = fill(fill_Xcoord, fillYcoord, light_green,'LineStyle','none');
            set(f,'facealpha',.5)
            
            legend([re,d_med, a_med, c_med],'RE', 'Decreasing', 'Anchor', 'Constant','location', 'southoutside')
            legend('boxoff')
        end
        title(titles(i))
        ax = gca; % current axes
        ax.FontSize = fs_prop;
        grid on
        grid minor
        %         if s==1
        %             ylim([-0.02, 0.01])
        %         elseif s==2
        %             ylim([-1, 0.4])
        %         elseif s==3
        %             ylim([-0.3, 0.1])
        %         end
    end
    if print_figs ==1
        figname = [this_code, '_', 'IRFs_' shocknames{s},'_rho',rho_val, '_psi_pi_', psi_pi_val]
        cd(figpath)
        export_fig(figname)
        cd(current_dir)
        close
    end
    
    % Plot gain conditional on shock
    figure
    set(gcf,'color','w'); % sets white background color
    set(gcf, 'Position', get(0, 'Screensize')); % sets the figure fullscreen
    mk = plot(med_k,'r','linewidth', 2); hold on
    mkc = plot(med_k_cusum,'color', saddle_brown,'linewidth', 2);
    fill_Xcoord2 = [1:T, fliplr(1:T)];
    fillYcoord2 = [lb_k', fliplr(ub_k')];
    f = fill(fill_Xcoord2, fillYcoord2, light_salmon,'LineStyle','none');
    set(f,'facealpha',.5)
    fillYcoord2 = [lb_k_cusum', fliplr(ub_k_cusum')];
    f = fill(fill_Xcoord2, fillYcoord2, light_brown,'LineStyle','none');
    set(f,'facealpha',.5)
    ax = gca; % current axes
    ax.FontSize = fs_prop;
    grid on
    grid minor
    legend([mk, mkc],'CEMP criterion',  'CUSUM criterion','location', 'southoutside')
    legend('boxoff')
    title('Inverse gain')
    if print_figs ==1
        figname = [this_code, '_', 'gain_',shocknames{s}, '_rho',rho_val, '_psi_pi_', psi_pi_val]
        cd(figpath)
        export_fig(figname)
        cd(current_dir)
        close
    end
    
    % IRF of the mean gain
    figure
    set(gcf,'color','w'); % sets white background color
    set(gcf, 'Position', get(0, 'Screensize')); % sets the figure fullscreen
    mk = plot(RIRk,'r','linewidth', 2); hold on
    mkc = plot(RIRk_cusum,'color', saddle_brown,'linewidth', 2);
    ax = gca; % current axes
    ax.FontSize = fs_prop;
    grid on
    grid minor
    legend([mk, mkc],'CEMP criterion',  'CUSUM criterion','location', 'southoutside')
    legend('boxoff')
    title('Inverse gain IRF (mean)')
    if print_figs ==1
        figname = [this_code, '_', 'gain_IRF_',shocknames{s}, '_rho',rho_val, '_psi_pi_', psi_pi_val]
        cd(figpath)
        export_fig(figname)
        cd(current_dir)
        close
    end
    
    figure
    set(gcf,'color','w'); % sets white background color
    set(gcf, 'Position', get(0, 'Screensize')); % sets the figure fullscreen
    for i=1:ny+1
        if i>ny
            subplot(1,ny+1,i)
            plot(zeros(1,h),'k--','linewidth', 2); hold on
%             fcst_RE = plot(RE_fcsts(1,:),'k','linewidth', 2);
            mf_c = plot(med_mf_c,'color', dark_green,'linewidth', 2);
        else
            subplot(1,ny+1,i)
            plot(zeros(1,h),'k--','linewidth', 2); hold on
            % Plot medians
            c_med = plot(medc(i,:),'color', dark_green,'linewidth', 2);
            % Plot CIs
            fill_Xcoord = [1:h, fliplr(1:h)];
            fillYcoord = [lbc(i,:), fliplr(ubc(i,:))];
            f = fill(fill_Xcoord, fillYcoord, light_green,'LineStyle','none');
            set(f,'facealpha',.5) 
        end
        title(titles(i))
        ax = gca; % current axes
        ax.FontSize = fs_prop;
        grid on
        grid minor
    end
    if print_figs ==1
        figname = [this_code, '_', 'IRFs_cgain_' shocknames{s},'_rho',rho_val, '_psi_pi_', psi_pi_val]
        cd(figpath)
        export_fig(figname)
        cd(current_dir)
        close
    end
    
end


disp(['(psi_x, psi_pi, rho)=   ', num2str([psi_x, psi_pi, rho])])
toc
