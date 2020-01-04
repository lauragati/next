% materials 6 - second try
% Goals:
% 1. simplify and unify learning codes and GIRs for learning
% 2. add interest rate smoothing w/o screwing up everything
% 3. implement IRFs conditional on when anchored or not
% 19 Oct 2019
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
skip_old_plots =0;

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
n = 3; % fortunately, this is the dimension of everything
P = eye(n).*[rho_r, rho_i, rho_u]';
SIG = eye(n).*[sig_r, sig_i, sig_u]';


% Sequence of innovations
rng(0)
e = randn(n,T);

% Solve RE model
[fyn, fxn, fypn, fxpn] = model_NK(param);
[gx,hx]=gx_hx_alt(fyn,fxn,fypn,fxpn);
[ny, nx] = size(gx);
[~, ~, Aa_LH, Ab_LH, As_LH] = matrices_A(param, hx);
[Aa, Ab, As] = matrices_A2(param, hx); % try with 3 Jan 2020 general method (MN) - you get the exact same matrices
% Aa_LH = Aa;
% Ab_LH = Ab;
% As_LH = As;

% Simulate models
% Use Ryan's code to simulate from the RE model
[x_RE, y_RE] = sim_model(gx,hx,SIG,T,burnin,e);

% LR learning (only constant, only for pi), decreasing gain
anal = 1; % take analytical LR exp
H=0;
[x_LR, y_LR] = sim_learn_LR_constant_pidrift(gx,hx,SIG,T,burnin,e, Aa_LH, Ab_LH, As_LH, param, setp, H, anal);

% LR learning with anchoring, learning the inflation drift only
[x_LR_anchor, y_LR_anchor, ~, ~, pibar, kd] = sim_learn_LR_anchoring_pidrift(gx,hx,SIG,T,burnin,e, Aa_LH, Ab_LH, As_LH, param, setp, H, anal);
kd_inv = 1./kd;

% LR learning (only constant, only for pi), constant gain ('perpetual learning')
[x_LR_perp, y_LR_perp] = sim_learn_LR_constant_pidrift_perpetual(gx,hx,SIG,T,burnin,e, Aa_LH, Ab_LH, As_LH, param, setp, H, anal);

% LR learning with anchoring, learning the inflation drift only, using
% alternative CUSUM criterion
[x_LR_anchor_cusum, y_LR_anchor_cusum, ~, ~, pibar_cusum, kd_cusum] = sim_learn_LR_anchoring_pidrift_cusum(gx,hx,SIG,T,burnin,e, Aa_LH, Ab_LH, As_LH, param, setp, H, anal);
kd_inv_cusum = 1./kd_cusum;

% General learning code for LH model - not doing Euler equation learning
% cause Preston won't allow us
gain = 2; % 1 = decreasing gain, 2 = endogenous gain, 3 = constant gain
criterion = 2; % 1= CEMP's, 2 = CUSUM
dt = 0; % when shock imposed. If zero, then no shock
x0 = zeros(n,1); % shock vector (delta)
constant_only = 1; % learning constant only
dgain = 1;
again = 2;
cgain = 3;
critCEMP=1;
critCUSUM=2;
free=1; % use versions of the code that are n-free (use hx instead of P)
not_free=0;
[x, y] = sim_learn(gx,hx,SIG,T,burnin,e, Aa_LH, Ab_LH, As_LH, param, constant_only, dgain, critCEMP);
% cool - I've verified that I'm getting the same things as before

% Gather observables
Y(:,:,1) =y_RE;
Y(:,:,2) =y_LR;
Y(:,:,3) =y_LR_anchor;
Y(:,:,4) =y_LR_perp;
Y(:,:,5) =y_LR_anchor_cusum;


if skip_old_plots==0
    %% Analysis plots
    % Observables
    figure
    set(gcf,'color','w'); % sets white background color
    set(gcf, 'Position', get(0, 'Screensize')); % sets the figure fullscreen
    
    titles = {'Inflation','Output gap','Int. rate'};
    l=0;
    for j=1:ny
        l = l+1;
        subplot(1,ny,l)
        plot(squeeze(Y(j,:,1)),'k','linewidth', 2); hold on
        plot(squeeze(Y(j,:,2)),'b','linewidth', 2)
        plot(squeeze(Y(j,:,3)),'r','linewidth', 2)
        plot(squeeze(Y(j,:,4)),'color', dark_green,'linewidth', 2)
        
        ax = gca; % current axes
        ax.FontSize = fs;
        grid on
        grid minor
        legend('RE', 'LR decreasing gain','LR anchor',  'LR constant gain')
        title(titles(l))
    end
    if print_figs ==1
        figname = [this_code, '_', 'observables_notfree']
        cd(figpath)
        export_fig(figname)
        cd(current_dir)
        close
    end
    
    % Gain
    figure
    set(gcf,'color','w'); % sets white background color
    set(gcf, 'Position', get(0, 'Screensize')); % sets the figure fullscreen
    
    subplot(1,2,1)
    plot(kd_inv(1,:),'r','linewidth', 2)
    ax = gca; % current axes
    ax.FontSize = fs;
    grid on
    grid minor
    title('Inverse gain')
    subplot(1,2,2)
    plot(pibar,'r','linewidth', 2)
    ax = gca; % current axes
    ax.FontSize = fs;
    grid on
    grid minor
    title('Inflation drift')
    if print_figs ==1
        figname = [this_code, '_', 'gain_drift_notfree']
        cd(figpath)
        export_fig(figname)
        cd(current_dir)
        close
    end
    
end
%% Generate IRFs for RE and anchoring
d = 1; % the innovation, delta
shocknames = {'natrate', 'monpol','costpush'};


% cycle thru the shocks of the model
for s=1:n
    x0 = zeros(1,n);
    x0(s) = d;
    h = 20; % horizon of IRF
    [IR, iry, irx]=ir(gx,hx,x0,h);
    iry = iry';
    
    % For learning, the approach I take is resimulate everything and for each
    % period t, expose the econ to this same shock. Then take an average.
    GIRd = zeros(n,h,T-h); % decreasing gain
    GIRa = zeros(n,h,T-h); % anchoring
    GIRc = zeros(n,h,T-h); % constant gain
    for t=1:T-h
        % 1. create alternative simulations, adding the impulse always at a new time t
        [x_LR_shockd, y_LR_shockd] = sim_learn_LR_constant_pidrift_shockd(gx,hx,SIG,T,burnin,e, Aa_LH, Ab_LH, As_LH, param, setp, H, anal, t, x0);
        [x_LR_anchor_shockd, y_LR_anchor_shockd] = sim_learn_LR_anchoring_pidrift_shockd(gx,hx,SIG,T,burnin,e, Aa_LH, Ab_LH, As_LH, param, setp, H, anal, t, x0);
        [x_LR_perp_shockd, y_LR_perp_shockd] = sim_learn_LR_constant_pidrift_perpetual_shockd(gx,hx,SIG,T,burnin,e, Aa_LH, Ab_LH, As_LH, param, setp, H, anal, t, x0);
        % Now let's shock our general learning code
        [xs, ys] = sim_learn(gx,hx,SIG,T,burnin,e, Aa_LH, Ab_LH, As_LH, param, constant_only, cgain, critCEMP, t, x0);
        %         y_LR_perp_shockd - y; % perfect - all of these are the same!
        % 2. take differences between this and the standard simulation
        GIRd(:,:,t) = y_LR_shockd(:,t:t+h-1) - y_LR(:,t:t+h-1);
        GIRa(:,:,t) = y_LR_anchor_shockd(:,t:t+h-1) - y_LR_anchor(:,t:t+h-1);
        GIRc(:,:,t) = y_LR_perp_shockd(:,t:t+h-1) - y_LR_perp(:,t:t+h-1);
    end
    
    % option 1: take simple averages
    RIRd = mean(GIRd,3);
    RIRa = mean(GIRa,3);
    RIRc = mean(GIRc,3);
    
    if skip_old_plots==0
        % Plot IRFs
        figure
        set(gcf,'color','w'); % sets white background color
        set(gcf, 'Position', get(0, 'Screensize')); % sets the figure fullscreen
        for i=1:n
            subplot(1,n,i)
            plot(zeros(1,h),'k--','linewidth', 2); hold on
            for j=2:h
                plot(squeeze(GIRd(i,:,j)),'color',light_sky_blue,'linewidth', 2)
                plot(squeeze(GIRa(i,:,j)),'color',light_salmon,'linewidth', 2)
                plot(squeeze(GIRc(i,:,j)),'color',light_green,'linewidth', 2)
            end
            re = plot(iry(i,:),'k','linewidth', 2);
            d_mean = plot(RIRd(i,:),'b','linewidth', 2);
            am_mean = plot(RIRa(i,:),'r','linewidth', 2);
            c_mean = plot(RIRc(i,:),'color', dark_green,'linewidth', 2);
            legend([re,d_mean, am_mean, c_mean],'RE', 'Decreasing gain', 'Anchor', 'Constant gain')
            title(titles(i))
            ax = gca; % current axes
            ax.FontSize = fs;
            grid on
            grid minor
            if s==1
                ylim([-0.015, 0.005])
            elseif s==3
                ylim([-0.2, 0.05])
            end
        end
        if print_figs ==1
            figname = [this_code, '_', 'IRFs_notfree_', shocknames{s}]
            cd(figpath)
            export_fig(figname)
            cd(current_dir)
            close
        end
    end
end
return
%% compare the two anchoring mechanisms
% close all
if skip_old_plots==0
    
    % Observables
    figure
    set(gcf,'color','w'); % sets white background color
    set(gcf, 'Position', get(0, 'Screensize')); % sets the figure fullscreen
    
    titles = {'Inflation','Output gap','Int. rate'};
    l=0;
    for j=1:ny
        l = l+1;
        subplot(1,ny,l)
        plot(squeeze(Y(j,:,3)),'r','linewidth', 2); hold on
        plot(squeeze(Y(j,:,5)),'color', saddle_brown,'linewidth', 2)
        ax = gca; % current axes
        ax.FontSize = fs;
        grid on
        grid minor
        legend('CEMP criterion',  'CUSUM criterion')
        title(titles(l))
    end
    if print_figs ==1
        figname = [this_code, '_', 'observables2_notfree']
        cd(figpath)
        export_fig(figname)
        cd(current_dir)
        close
    end
    
    % Gain
    figure
    set(gcf,'color','w'); % sets white background color
    set(gcf, 'Position', get(0, 'Screensize')); % sets the figure fullscreen
    
    subplot(1,2,1)
    plot(kd_inv(1,:),'r','linewidth', 2); hold on
    plot(kd_inv_cusum(1,:),'color', saddle_brown,'linewidth', 2)
    ax = gca; % current axes
    ax.FontSize = fs;
    grid on
    grid minor
    legend('CEMP criterion',  'CUSUM criterion')
    title('Inverse gain')
    subplot(1,2,2)
    plot(pibar,'r','linewidth', 2); hold on
    plot(pibar_cusum,'color', saddle_brown,'linewidth', 2)
    ax = gca; % current axes
    ax.FontSize = fs;
    grid on
    grid minor
    legend('CEMP criterion',  'CUSUM criterion')
    title('Inflation drift')
    if print_figs ==1
        figname = [this_code, '_', 'gain_drift2_notfree']
        cd(figpath)
        export_fig(figname)
        cd(current_dir)
        close
    end
    
end

close all


%% New approach:  implement interest rate smoothing AFTER you've done the rest
% This way I'm hoping to preserve things intact
% close all

% introduce the new parameter and adapt shit
rho = param.rho;
% introduce adaptive names depending on the value of rho
rho_val_raw = num2str(rho);
rho_val = replace(rho_val_raw,'.','_');
psi_pi_val_raw = num2str(psi_pi);
psi_pi_val = replace(psi_pi_val_raw,'.','_');

n = 4; % now n becomes 4
% P = eye(n).*[rho_r, rho_i, rho_u, 0]';
SIG = eye(n).*[sig_r, sig_i, sig_u, 0]';


% Sequence of innovations
% adopt the old one - somehow, rng(0)-ing doesn't get you the same thing
% because e_i is a different size
e_i = [e; zeros(1,T)];
% gzactly - note that the below reproduces e
% rng(0)
% e2 = randn(n-1,T);

% Solve RE model
[fyn, fxn, fypn, fxpn] = model_NK_intrate_smoothing(param);
[gx,hx]=gx_hx_alt(fyn,fxn,fypn,fxpn) % cool is the same so far
[ny, nx] = size(gx);
[~, ~, Aa_LH_i, Ab_LH_i, As_LH_i] = matrices_A_intrate_smoothing(param, setp, hx); % perfect - Aa, Ab are the same as before, As has extra column of zeros, otherwise identical!

% Simulate models
% Use Ryan's code to simulate from the RE model
[x_RE_i, y_RE_i] = sim_model(gx,hx,SIG,T,burnin,e_i);


% Now simulate the learning models using the general code
[x_d, y_d] = sim_learn(gx,hx,SIG,T,burnin,e_i, Aa_LH_i, Ab_LH_i, As_LH_i, param, setp, H, anal, constant_only, dgain, critCEMP, free);
[x_a, y_a, ~, ~,pibar_a, k_a] = sim_learn(gx,hx,SIG,T,burnin,e_i, Aa_LH_i, Ab_LH_i, As_LH_i, param, setp, H, anal, constant_only, again, critCEMP, free);
[x_a_cusum, y_a_cusum, ~, ~,pibar_a_cusum, k_cusum] = sim_learn(gx,hx,SIG,T,burnin,e_i, Aa_LH_i, Ab_LH_i, As_LH_i, param, setp, H, anal, constant_only, again, critCUSUM, free);
[x_c, y_c] = sim_learn(gx,hx,SIG,T,burnin,e_i, Aa_LH_i, Ab_LH_i, As_LH_i, param, setp, H, anal, constant_only, cgain, critCEMP, free);
ka_inv = 1./k_a;
k_cusum_inv = 1./k_cusum;


% Gather observables
Z(:,:,1) =y_RE_i;
Z(:,:,2) =y_d;
Z(:,:,3) =y_a;
Z(:,:,4) =y_c;
Z(:,:,5) =y_a_cusum;

%% Analysis plots interest rate smoothing
if skip_old_plots==0
    
    % Observables
    figure
    set(gcf,'color','w'); % sets white background color
    set(gcf, 'Position', get(0, 'Screensize')); % sets the figure fullscreen
    
    titles = {'Inflation','Output gap','Int. rate'};
    l=0;
    for j=1:ny
        l = l+1;
        subplot(1,ny,l)
        plot(squeeze(Z(j,:,1)),'k','linewidth', 2); hold on
        plot(squeeze(Z(j,:,2)),'b','linewidth', 2)
        plot(squeeze(Z(j,:,3)),'r','linewidth', 2)
        plot(squeeze(Z(j,:,4)),'color', dark_green,'linewidth', 2)
        
        ax = gca; % current axes
        ax.FontSize = fs_prop;
        grid on
        grid minor
        legend('RE', 'LH decreasing gain','LH anchor',  'LH constant gain', 'location', 'southoutside')
        title(titles(l))
    end
    if print_figs ==1
        figname = [this_code, '_', 'observables_intrate_smoothing_rho',rho_val, '_psi_pi_', psi_pi_val]
        cd(figpath)
        export_fig(figname)
        cd(current_dir)
        close
    end
    
    % Gain
    figure
    set(gcf,'color','w'); % sets white background color
    set(gcf, 'Position', get(0, 'Screensize')); % sets the figure fullscreen
    
    subplot(1,2,1)
    plot(ka_inv(1,:),'r','linewidth', 2)
    ax = gca; % current axes
    ax.FontSize = fs;
    grid on
    grid minor
    title('Inverse gain')
    subplot(1,2,2)
    plot(pibar_a,'r','linewidth', 2)
    ax = gca; % current axes
    ax.FontSize = fs_prop;
    grid on
    grid minor
    title('Inflation drift')
    if print_figs ==1
        figname = [this_code, '_', 'gain_drift_intrate_smoothing_rho',rho_val, '_psi_pi_', psi_pi_val]
        cd(figpath)
        export_fig(figname)
        cd(current_dir)
        close
    end
    
    % Ok finally! The graph is the same, and thu w/ interest rate smoothing, if
    % rho = 0, I obtain the exact same dynamics. Cool!
end
%% GIRs for interest rate smoothing
d = 1; % the innovation, delta
shocknames = {'natrate', 'monpol','costpush'};


% cycle thru the shocks of the model
for s=1:n-1
    x0 = zeros(1,n);
    x0(s) = d;
    h = 10; % horizon of IRF
    [IR, iry, irx]=ir(gx,hx,x0,h);
    iry = iry';
    
    % For learning, the approach I take is resimulate everything and for each
    % period t, expose the econ to this same shock. Then take an average.
    GIRd = zeros(ny,h,T-h); % decreasing gain
    GIRa = zeros(ny,h,T-h); % anchoring
    GIRc = zeros(ny,h,T-h); % constant gain
    pibars_a = zeros(T,1,T-h);
    pibars_a_cusum = zeros(T,1,T-h);
    ks_a = zeros(T,1,T-h);
    ks_a_cusum = zeros(T,1,T-h);
    for t=1:T-h
        % 1. create alternative simulations, adding the impulse always at a new time t
        % Now let's shock our general learning code
        [xs_d, ys_d] = sim_learn(gx,hx,SIG,T,burnin,e_i, Aa_LH_i, Ab_LH_i, As_LH_i, param, setp, H, anal, constant_only, dgain, critCEMP,free, t, x0);
        [xs_a, ys_a, ~, ~,pibars_a(:,:,t), ks_a(:,:,t)] = sim_learn(gx,hx,SIG,T,burnin,e_i, Aa_LH_i, Ab_LH_i, As_LH_i, param, setp, H, anal, constant_only, again, critCEMP,free, t, x0);
        [xs_c, ys_c] = sim_learn(gx,hx,SIG,T,burnin,e_i, Aa_LH_i, Ab_LH_i, As_LH_i, param, setp, H, anal, constant_only, cgain, critCEMP,free, t, x0);
        [xs_a_cusum, ys_a_cusum, ~, ~,pibars_a_cusum(:,:,t), ks_a_cusum(:,:,t)] = sim_learn(gx,hx,SIG,T,burnin,e_i, Aa_LH_i, Ab_LH_i, As_LH_i, param, setp, H, anal, constant_only, again, critCUSUM,free, t, x0);
        
        % 2. take differences between this and the standard simulation
        GIRd(:,:,t) = ys_d(:,t:t+h-1) - y_d(:,t:t+h-1);
        GIRa(:,:,t) = ys_a(:,t:t+h-1) - y_a(:,t:t+h-1);
        GIRc(:,:,t) = ys_c(:,t:t+h-1) - y_c(:,t:t+h-1);
    end
    
    % option 1: take simple averages
    RIRd = mean(GIRd,3);
    RIRa = mean(GIRa,3);
    RIRc = mean(GIRc,3);
    
    % take averages of gains and drifts too
    pibars_a_mean = mean(pibars_a,3);
    pibars_a_cusum_mean = mean(pibars_a_cusum,3);
    ks_a_mean = mean(ks_a,3);
    ks_a_cusum_mean = mean(ks_a_cusum,3);
    inv_ks_a_mean = 1./ks_a_mean;
    inv_ks_a_cusum_mean = 1./ks_a_cusum_mean;
    
    if skip_old_plots==0
        % Plot IRFs
        figure
        set(gcf,'color','w'); % sets white background color
        set(gcf, 'Position', get(0, 'Screensize')); % sets the figure fullscreen
        for i=1:n-1
            subplot(1,n-1,i)
            plot(zeros(1,h),'k--','linewidth', 2); hold on
            %         for j=2:h
            %             plot(squeeze(GIRd(i,:,j)),'color',light_sky_blue,'linewidth', 2)
            %             plot(squeeze(GIRa(i,:,j)),'color',light_salmon,'linewidth', 2)
            %             plot(squeeze(GIRc(i,:,j)),'color',light_green,'linewidth', 2)
            %         end
            re = plot(iry(i,:),'k','linewidth', 2);
            d_mean = plot(RIRd(i,:),'b','linewidth', 2);
            am_mean = plot(RIRa(i,:),'r','linewidth', 2);
            c_mean = plot(RIRc(i,:),'color', dark_green,'linewidth', 2);
            legend([re,d_mean, am_mean, c_mean],'RE', 'Decreasing gain', 'Anchor', 'Constant gain','location', 'southoutside')
            title(titles(i))
            ax = gca; % current axes
            ax.FontSize = fs_prop;
            grid on
            grid minor
            if s==1
                ylim([-0.02, 0.01])
            elseif s==2
                ylim([-1, 0.4])
            elseif s==3
                ylim([-0.3, 0.1])
            end
        end
        if print_figs ==1
            figname = [this_code, '_', 'IRFs_intrate_smoothing_', shocknames{s},'_rho',rho_val, '_psi_pi_', psi_pi_val]
            cd(figpath)
            export_fig(figname)
            cd(current_dir)
            close
        end
        
        % Plot gain and drift conditional on shock
        figure
        set(gcf,'color','w'); % sets white background color
        set(gcf, 'Position', get(0, 'Screensize')); % sets the figure fullscreen
        
        subplot(1,2,1)
        plot(inv_ks_a_mean,'r','linewidth', 2); hold on
        plot(inv_ks_a_cusum_mean,'color', saddle_brown,'linewidth', 2)
        ax = gca; % current axes
        ax.FontSize = fs_prop;
        grid on
        grid minor
        legend('CEMP criterion',  'CUSUM criterion','location', 'southoutside')
        title('Inverse gain')
        subplot(1,2,2)
        plot(pibars_a_mean,'r','linewidth', 2); hold on
        plot(pibars_a_cusum_mean,'color', saddle_brown,'linewidth', 2)
        ax = gca; % current axes
        ax.FontSize = fs_prop;
        grid on
        grid minor
        legend('CEMP criterion',  'CUSUM criterion','location', 'southoutside')
        title('Inflation drift')
        if print_figs ==1
            figname = [this_code, '_', 'gain_drift_',shocknames{s}, '_rho',rho_val, '_psi_pi_', psi_pi_val]
            cd(figpath)
            export_fig(figname)
            cd(current_dir)
            close
        end
    end
    
end

%% compare the two anchoring mechanisms with interest rate smoothing
% close all
if skip_old_plots==0
    
    % Observables
    figure
    set(gcf,'color','w'); % sets white background color
    set(gcf, 'Position', get(0, 'Screensize')); % sets the figure fullscreen
    
    titles = {'Inflation','Output gap','Int. rate'};
    l=0;
    for j=1:ny
        l = l+1;
        subplot(1,ny,l)
        plot(squeeze(Z(j,:,3)),'r','linewidth', 2); hold on
        plot(squeeze(Z(j,:,5)),'color', saddle_brown,'linewidth', 2)
        ax = gca; % current axes
        ax.FontSize = fs;
        grid on
        grid minor
        legend('CEMP criterion',  'CUSUM criterion')
        title(titles(l))
    end
    if print_figs ==1
        figname = [this_code, '_', 'observables_cusum_intrate_smoothing_rho',rho_val, '_psi_pi_', psi_pi_val]
        cd(figpath)
        export_fig(figname)
        cd(current_dir)
        close
    end
    
    % Gain
    figure
    set(gcf,'color','w'); % sets white background color
    set(gcf, 'Position', get(0, 'Screensize')); % sets the figure fullscreen
    
    subplot(1,2,1)
    plot(ka_inv(1,:),'r','linewidth', 2); hold on
    plot(k_cusum_inv(1,:),'color', saddle_brown,'linewidth', 2)
    ax = gca; % current axes
    ax.FontSize = fs_prop;
    grid on
    grid minor
    legend('CEMP criterion',  'CUSUM criterion','location', 'southoutside')
    title('Inverse gain')
    subplot(1,2,2)
    plot(pibar_a,'r','linewidth', 2); hold on
    plot(pibar_a_cusum,'color', saddle_brown,'linewidth', 2)
    ax = gca; % current axes
    ax.FontSize = fs_prop;
    grid on
    grid minor
    legend('CEMP criterion',  'CUSUM criterion','location', 'southoutside')
    title('Inflation drift')
    if print_figs ==1
        figname = [this_code, '_', 'gain_drift_cusum_intrate_smoothing_rho',rho_val, '_psi_pi_', psi_pi_val]
        cd(figpath)
        export_fig(figname)
        cd(current_dir)
        close
    end
    
    disp(['(psi_x, psi_pi, rho)=   ', num2str([psi_x, psi_pi, rho])])
end
%% Generating IRFs differentiating ones that are anchored / unanchored when shock hits
% let 1 denote anchored
% let 0 (or 2) denote unanchored
anchored = nan(1,T-h);
close all

d = 1; % the innovation, delta
shocknames = {'natrate', 'monpol','costpush'};
titles = {'Inflation','Output gap','Int. rate'};
% cycle thru the shocks of the model
for s=1:n-1
    x0 = zeros(1,n);
    x0(s) = d;
    h = 10; % horizon of IRF
    [IR, iry, irx]=ir(gx,hx,x0,h);
    iry = iry';
    
    % For learning, the approach I take is resimulate everything and for each
    % period t, expose the econ to this same shock. Then take an average.
    GIR = zeros(ny,h,T-h); % third index is when the shock is imposed
    pibars = zeros(T,1,T-h);
    ks = zeros(T,1,T-h);
    for t=1:T-h
        % 1. create alternative simulations, adding the impulse always at a new time t
        % Now let's shock our general learning code
        [xs, ys, ~, ~,pibars(:,:,t), ks(:,:,t), anchored(t)] = sim_learn(gx,hx,SIG,T,burnin,e_i, Aa_LH_i, Ab_LH_i, As_LH_i, param, setp, H, anal, constant_only, again, critCEMP,free, t, x0);
        
        % 2. take differences between this and the standard simulation
        GIR(:,:,t) = ys(:,t:t+h-1) - y_a(:,t:t+h-1);
    end
    
    % get indeces for anchored or unanchored at time shock hits
    anch = find(anchored==1);
    unanch = find(anchored==0);
    % option 1: take simple averages
    RIR1 = mean(GIR(:,:,anch),3);
    RIR2 = mean(GIR(:,:,unanch),3);
    % take averages of gains and drifts too
    pibar_mean1 = mean(pibars(:,:,anch),3);
    pibar_mean2 = mean(pibars(:,:,unanch),3);
    ks_mean1 = mean( ks(:,:,anch),3);
    ks_mean2 = mean( ks(:,:,unanch),3);
    inv_ks_mean1 = 1./ks_mean1;
    inv_ks_mean2 = 1./ks_mean2;
    
    
    % option 2: sort and take percentile bands
    [lb1, med1, ub1] = confi_bands(GIR(:,:,anch),0.1);
    [lb2, med2, ub2] = confi_bands(GIR(:,:,unanch),0.1);
    
    [lb_k1, med_k1, ub_k1] = confi_bands(1./ks(:,:,anch),0.1);
    [lb_k2, med_k2, ub_k2] = confi_bands(1./ks(:,:,unanch),0.1);
    
    [lb_pib1, med_pib1, ub_pib1] = confi_bands(pibars(:,:,anch),0.1);
    [lb_pib2, med_pib2, ub_pib2] = confi_bands(pibars(:,:,unanch),0.1);
    
    % Plot IRFs
    
    % IRF of Option 1
    skip_this=1;
    if skip_this==0
        figure
        set(gcf,'color','w'); % sets white background color
        set(gcf, 'Position', get(0, 'Screensize')); % sets the figure fullscreen
        for i=1:n-1
            subplot(1,n-1,i)
            plot(zeros(1,h),'k--','linewidth', 2); hold on
            %         for j=2:h
            %             plot(squeeze(GIRd(i,:,j)),'color',light_sky_blue,'linewidth', 2)
            %             plot(squeeze(GIRa(i,:,j)),'color',light_salmon,'linewidth', 2)
            %             plot(squeeze(GIRc(i,:,j)),'color',light_green,'linewidth', 2)
            %         end
            re = plot(iry(i,:),'k','linewidth', 2);
            ancho = plot(RIR1(i,:),'b','linewidth', 2);
            unancho = plot(RIR2(i,:),'r--','linewidth', 2);
            legend([re, ancho, unancho], 'RE', 'Anchored', 'Unanchored','location', 'southoutside')
            title(titles(i))
            ax = gca; % current axes
            ax.FontSize = fs_prop;
            grid on
            grid minor
            if s==1
                ylim([-0.02, 0.01])
            elseif s==2
                ylim([-1, 0.4])
            elseif s==3
                ylim([-0.3, 0.1])
            end
        end
        if print_figs ==1
            figname = [this_code, '_', 'IRFs_intrate_smoothing_', 'cond_anch_', shocknames{s},'_rho',rho_val, '_psi_pi_', psi_pi_val]
            cd(figpath)
            export_fig(figname)
            cd(current_dir)
            close
        end
        close
    end
    % IRF of Option 2
    figure
    set(gcf,'color','w'); % sets white background color
    set(gcf, 'Position', get(0, 'Screensize')); % sets the figure fullscreen
    for i=1:n-1
        subplot(1,n-1,i)
        plot(zeros(1,h),'k--','linewidth', 2); hold on
        re = plot(iry(i,:),'k','linewidth', 2);
        ancho = plot(med1(i,:),'b','linewidth', 2);
        fill_Xcoord = [1:h, fliplr(1:h)];
        fillYcoord = [lb1(i,:), fliplr(ub1(i,:))];
        f = fill(fill_Xcoord, fillYcoord, light_sky_blue,'LineStyle','none');
        set(f,'facealpha',.5)
        
        unancho = plot(med2(i,:),'r--','linewidth', 2);
        fillYcoord = [lb2(i,:), fliplr(ub2(i,:))];
        f = fill(fill_Xcoord, fillYcoord, light_salmon,'LineStyle','none');
        set(f,'facealpha',.5)
        uistack(ancho, 'top') % move f to the top layer of the figure
        uistack(unancho, 'top') % move f to the top layer of the figure
        uistack(re, 'top') % move f to the top layer of the figure
        legend([re, ancho, unancho], 'RE', 'Anchored', 'Unanchored','location', 'southoutside')
        title(titles(i))
        ax = gca; % current axes
        ax.FontSize = fs_prop;
        grid on
        grid minor
        if s==1
            ylim([-0.02, 0.01])
        elseif s==2
            ylim([-1, 0.4])
        elseif s==3
            ylim([-0.3, 0.1])
        end
    end
    if print_figs ==1
        figname = [this_code, '_', 'IRFs_intrate_smoothing_', 'cond_anch_','SORTED_', shocknames{s},'_rho',rho_val, '_psi_pi_', psi_pi_val]
        cd(figpath)
        export_fig(figname)
        cd(current_dir)
        close
    end
    
    
    % Plot gain and drift conditional on shock
    figure
    set(gcf,'color','w'); % sets white background color
    set(gcf, 'Position', get(0, 'Screensize')); % sets the figure fullscreen
    plot(inv_ks_mean1,'b','linewidth', 2); hold on
    plot(inv_ks_mean2,'r--','linewidth', 2)
    fill_Xcoord = [1:T, fliplr(1:T)];
    fillYcoord = [lb_k1', fliplr(ub_k1)'];
    f = fill(fill_Xcoord, fillYcoord, light_sky_blue,'LineStyle','none');
    set(f,'facealpha',.5)
    fillYcoord = [lb_k2', fliplr(ub_k2)'];
    f = fill(fill_Xcoord, fillYcoord, light_salmon,'LineStyle','none');
    set(f,'facealpha',.5)
    ax = gca; % current axes
    ax.FontSize = fs_prop;
    grid on
    grid minor
    legend('Anchored',  'Unanchored','location', 'southoutside')
    title('Inverse gain')
    if print_figs ==1
        figname = [this_code, '_', 'cond_anch_','gain_',shocknames{s}, '_rho',rho_val, '_psi_pi_', psi_pi_val]
        cd(figpath)
        export_fig(figname)
        cd(current_dir)
        close
    end
    skip_this=1;
    if skip_this==0
        figure
        set(gcf,'color','w'); % sets white background color
        set(gcf, 'Position', get(0, 'Screensize')); % sets the figure fullscreen
        plot(pibar_mean1,'b','linewidth', 2); hold on
        plot(pibar_mean2,'r--','linewidth', 2)
        ax = gca; % current axes
        ax.FontSize = fs_prop;
        grid on
        grid minor
        legend('Anchored',  'Unanchored','location', 'southoutside')
        title('Inflation drift')
        if print_figs ==1
            figname = [this_code, '_', 'cond_anch_','drift_',shocknames{s}, '_rho',rho_val, '_psi_pi_', psi_pi_val]
            cd(figpath)
            export_fig(figname)
            cd(current_dir)
            close
        end
        close
    end
end

disp(['(psi_x, psi_pi, rho)=   ', num2str([psi_x, psi_pi, rho])])
disp('Done.')
toc