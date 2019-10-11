% materials 5
% work after DW
% 7 Oct 2019
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
lw=2; % linewidth
fs_pres = 80;
lw_pres = 6;
fs_prop = 40;
lw_prop = 4;

%% Simulation
T = 40
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
[Ap_RE, As_RE, Aa_LR, Ab_LR, As_LR, B1, B2] = matrices_A(param, setp);

% Simulate the three models
% Use Ryan's code to simulate from the RE model
[x_RE, y_RE] = sim_model(gx,hx,SIG,T,burnin,e);

% LR learning (only constant, only for pi) 
anal = 1; % take analytical LR exp
H=0;
[x_LR, y_LR] = sim_learn_LR_constant_pidrift(gx,hx,SIG,T,burnin,e, Aa_LR, Ab_LR, As_LR, param, setp, H, anal);

% LR learning with anchoring, learning the inflation drift only
anal = 1; % take analytical LR exp
[x_LR_anchor, y_LR_anchor, ~, ~, pibar, kd] = sim_learn_LR_anchoring_pidrift(gx,hx,SIG,T,burnin,e, Aa_LR, Ab_LR, As_LR, param, setp, H, anal);
kd_inv = 1./kd;


% Gather observables
Y(:,:,1) =y_RE; % let's call em V for a change
Y(:,:,2) =y_LR;
Y(:,:,3) =y_LR_anchor;

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
    plot(squeeze(Y(j,:,1)),'linewidth', 2); hold on
    plot(squeeze(Y(j,:,2)),'linewidth', 2)
    plot(squeeze(Y(j,:,3)),'k','linewidth', 2)
    ax = gca; % current axes
    ax.FontSize = fs;
    grid on
    grid minor
    legend('RE', 'LR','LR anchor')
    title(titles(l))
end
if print_figs ==1
    figname = ['materials5_observables1']
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
plot(kd_inv(1,:),'linewidth', 2)
ax = gca; % current axes
ax.FontSize = fs;
grid on
grid minor
title('Inverse gain')
subplot(1,2,2)
plot(pibar,'linewidth', 2)
ax = gca; % current axes
ax.FontSize = fs;
grid on
grid minor
title('Inflation drift')
if print_figs ==1
    figname = ['materials5_gain_drift']
    cd(figpath)
    export_fig(figname)
    cd(current_dir)
    close
end


%% Generate IRFs for RE and anchoring
x0 = [0 1 0]'; % shock to monpol
h = 10; % horizon of IRF
[IR, iry, irx]=ir(gx,hx,x0,h);

% For learning, the approach I take is resimulate everything and for each
% period t, expose the econ to this same shock. Then take an average.

GIR = zeros(n,h,T-h);
for t=1:T-h
    % 1. create alternative simulations, adding the impulse always at a new
    % time t
    [x_LR_anchor_shockd, y_LR_anchor_shockd] = sim_learn_LR_anchoring_pidrift_shockd(gx,hx,SIG,T,burnin,e, Aa_LR, Ab_LR, As_LR, param, setp, H, anal, t);
    % 2. take differences between this and the standard simulation
    GIR(:,:,t) = y_LR_anchor(:,t:t+h-1) - y_LR_anchor_shockd(:,t:t+h-1);
end

% option 1: take simple averages
RIR1 = mean(GIR,3);
