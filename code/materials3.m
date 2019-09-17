% materials 3
% Simulating special cases of the mixed CEMP-Preston model
% 12 Sep 2019
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

% Variable stuff ---
print_figs    = 0;
do_old_plots  = 0;
if print_figs ==1
    output_table  = 1;
else
    output_table =0;
end

fs=20; % fontsize


%% Simulation
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

tic
T = 2000000;
burnin = 0;
% learning convergence tolerance
tol= 1e-4;

% Sequence of shocks
rng(0)
e = randn(n,T);
s = zeros(n,T);
s(:,1) = e(:,1);
for t=2:T
    s(:,t) = P*s(:,t-1) + SIG*e(:,t);
end

% Solve RE model
[fyn, fxn, fypn, fxpn] = model_NK(param);
[gx,hx]=gx_hx_alt(fyn,fxn,fypn,fxpn);
[ny, nx] = size(gx);
[Ap_RE, As_RE, Aa_LR, Ab_LR, As_LR, B1, B2] = matrices_A(param, setp);
%
top1 = -[fyn fxpn]\(fypn);
Ap_RE - top1(1:3,:); % yep they equal
top2 = -[fyn fxpn]\(fxn);
As_RE - top2(1:3,:); % yep they equal

H = 10000;
% phi_rand = rand(3,4);
% [fa_anal, fb_anal] = fafb_anal(param, setp, phi_rand, s(:,1));
% [fa_trunc, fb_trunc] = fafb_trunc(param, setp, phi_rand, s(:,1), H);
% fa_anal - fa_trunc % these are the same
% fb_anal - fb_trunc % these are the same (to -12) only if H>10000.


% Simulate the three models
% Use Ryan's code to simulate from the RE model
[x_RE, y_RE] = sim_model(gx,hx,SIG,T,burnin,e); % NOTE: need to input innovations here such that x = s!
% Note:
% y = z
% x = s

% Euler equation learning (learning both slope and constant)
[x_EE, y_EE] = sim_learn_EE(gx,hx,fxpn,fxn,fypn,fyn,SIG,T,burnin,e);

% Note:
% x_RE = x_EE = s
% b/c 1) agents aren't learning about states 2) in NK
% model, the only states are exogenous states

% a check (learning both slope and constant)
[x_EE2, y_EE2] = sim_learn_EE_check(gx,hx,fxpn,fxn,fypn,fyn,SIG,T,burnin,e, Ap_RE, As_RE, param, setp);
% x_EE = x_EE2
y_EE - y_EE2; % yep they equal too


% LR learning (learning both slope and constant)
anal = 1; % take analytical LR exp
[x_LR, y_LR, ~, diff] = sim_learn_LR(gx,hx,SIG,T,burnin,e, Aa_LR, Ab_LR, As_LR, param, setp, H, anal);

% Gather observables
Y(:,:,1) =y_RE;
Y(:,:,2) =y_EE;
Y(:,:,3) =y_LR;

skip_this=1;
if skip_this==0
    % Now do an exercise in which horizons of expectations are compared
    Z(:,:,1) = y_LR;
    anal = 0; % take analytical LR exp
    horizons = [10000, 100, 1];
    for i=1:3
        H = horizons(i);
        tic
        [~, Z(:,:,i+1)] = sim_learn_LR(gx,hx,SIG,T,burnin,e, Aa_LR, Ab_LR, As_LR, param, setp, H, anal);
        toc
    end
    Z(:,:,end+1) = y_EE;
end
%% Plots
% Observables
figure
set(gcf,'color','w'); % sets white background color
set(gcf, 'Position', get(0, 'Screensize')); % sets the figure fullscreen

titles = {'Inflation','Output gap','Int. rate'};
k=0;
for j=1:ny
    k = k+1;
    subplot(1,ny,k)
    plot(squeeze(Y(j,:,1)),'linewidth', 2); hold on
    plot(squeeze(Y(j,:,2)),'linewidth', 2)
    plot(squeeze(Y(j,:,3)),'linewidth', 2)
    ax = gca; % current axes
    ax.FontSize = fs;
    grid on
    grid minor
    legend('RE', 'EE','LR')
    title(titles(k))
end
if print_figs ==1
    figname = ['materials3_observables1']
    cd(figpath)
    export_fig(figname)
    cd(current_dir)
    close
end
% close

if skip_this==0
    % Comparing horizons
    figure
    set(gcf,'color','w'); % sets white background color
    set(gcf, 'Position', get(0, 'Screensize')); % sets the figure fullscreen
    
    titles = {'Inflation','Output gap','Int. rate'};
    k=0;
    for j=1:ny
        k = k+1;
        subplot(2,ny,k)
        plot(squeeze(Z(j,:,1)),'linewidth', 2); hold on
        plot(squeeze(Z(j,:,2)),'linewidth', 2, 'linestyle', '--')
        plot(squeeze(Z(j,:,3)),'linewidth', 2)
        ax = gca; % current axes
        ax.FontSize = fs;
        grid on
        grid minor
        legend('LR analytical', 'H = 10000','H = 100', 'location', 'best')
        title(titles(k))
        subplot(2,ny,k+3)
        %     plot(squeeze(Z(j,:,3)),'linewidth', 2)
        plot(squeeze(Z(j,:,4)),'linewidth', 2); hold on
        plot(squeeze(Z(j,:,5)),'linewidth', 2)
        ax = gca; % current axes
        ax.FontSize = fs;
        grid on
        grid minor
        legend('H = 1', 'EE', 'location', 'best')
        title(titles(k))
    end
    if print_figs ==1
        figname = ['materials3_horizons']
        cd(figpath)
        export_fig(figname)
        cd(current_dir)
        close
    end
end

% Max abs differences in phi, LR learning
figure
set(gcf,'color','w'); % sets white background color
set(gcf, 'Position', get(0, 'Screensize')); % sets the figure fullscreen
plot(diff,'linewidth', 2)
ylim([0, 0.05])
ax = gca; % current axes
ax.FontSize = fs;
grid on
grid minor
title('Maximum absolute differences in \phi, LR model')
if print_figs ==1
    figname = ['materials3_diffs_LR']
    cd(figpath)
    export_fig(figname)
    cd(current_dir)
    close
end

% RE and LR models only, last 100 periods
figure
set(gcf,'color','w'); % sets white background color
set(gcf, 'Position', get(0, 'Screensize')); % sets the figure fullscreen
plot(diff,'linewidth', 2)
k=0;
for j=1:ny
    k = k+1;
    subplot(1,ny,k)
    plot(squeeze(Y(j,end-100:end,1)),'linewidth', 2); hold on
    plot(squeeze(Y(j,end-100:end,3)),'linewidth', 2)
    ax = gca; % current axes
    ax.FontSize = fs;
    grid on
    grid minor
    legend('RE','LR', 'location', 'best')
    title(titles(k))
end
if print_figs ==1
    figname = ['materials3_observables_RE_LR_last100']
    cd(figpath)
    export_fig(figname)
    cd(current_dir)
    close
end

toc