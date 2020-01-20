% materials 13
% Goals: (still!)
% Find the changes in the model that make it fit data
% 1.) Changes to the policy rule
% 2.) Changes to expectation formation
% 16 Jan 2020
clearvars
close all
clc

% Add all the relevant paths and grab the codename
this_code = mfilename;
[current_dir, basepath, BC_researchpath,toolpath,export_figpath,figpath,tablepath,datapath] = add_paths;

% Variable stuff ---
print_figs        = 0;
stop_before_plots = 0;
skip_old_plots    = 0;
output_table = print_figs;

skip_old_stuff = 1;


%% Do IRFs for extensions - same code as before
command_IRFs_many_learning

%% Trying to do an analog of the Kalman gain

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

% Take the very baseline model
[fyn, fxn, fypn, fxpn] = model_NK(param);
[gx,hx]=gx_hx_alt(fyn,fxn,fypn,fxpn);
[ny, nx] = size(gx);
SIG = eye(nx).*[sig_r, sig_i, sig_u]';
eta = SIG;
[~,sigx] = mom(gx,hx,eta*eta');
kapp_star = sigx*gx'*(gx*sigx*gx')^(-1)

%% A simple LOM of regime r
% P = [p11, p12; p21, p22];
% Follow Davig and Leeper's transition probabilities:
% and their notation of 1 being the active regime (psi_pi = 2.19)
% and 2 the passive regime (psi_pi = 0.89).
p11 = 0.95; %0.95
p21 = 1-p11;
p22 = 0.93; % 0.93
p12 = 1-p22;

r = nan(1,10);
r(1) = 1; % initialize in the active regime
for t=2:20
    if r(t-1)==1 % active regime
    r(t) = randsample([1,2],1,'true',[p11, p21]);
    elseif r(t-1)==2 % passive regime
    r(t) = randsample([1,2],1,'true',[p12, p22]);
    end
end
r % it works, tho a little messy
