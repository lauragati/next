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

T = 500;
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
R0 = eye(n);
phi0 = zeros(n,1);


[Ap_RE, As_RE, Aa_LR, Ab_LR, As_LR, B1, B2] = matrices_A(param, setp);

H = 100;
[fa_anal, fb_anal] = fafb_anal(param, setp, phi0, s(:,1));
[fa_trunc, fb_trunc] = fafb_trunc(param, setp, phi0, s(:,1), H); % they are equal, cool


% Simulate the three models
% Use Ryan's code to simulate from the RE model
[x1, y1] = sim_model(gx,hx,SIG,T,burnin,e); % NOTE: need to input innovations here such that x1 = s!

return
% Euler equation learning
tic
[z_EE, conv_EE] = sim_learn_EE(param,gx,R0,phi0,Ap_RE,As_RE,s,tol);
toc

% LR learning
tic
[zz] = sim_learn_LR(param,R0,phi0,s,tol); 
toc



