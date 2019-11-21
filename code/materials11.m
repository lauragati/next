% materials 11
% Goals:
% 1. find a reasonable value of the gain gbar
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

skip_old_stuff = 1;

%% Simulation
tic
T = 400 % 400
burnin = 0;

% Cross-section
N = 500; %500 size of cross-section
eN = randn(ne,T,N); % gen all the N sequences of shocks at once.

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


% Model with interest rate smoothing
[fyn, fxn, fypn, fxpn] = model_NK_intrate_smoothing(param);
[gx,hx]=gx_hx_alt(fyn,fxn,fypn,fxpn);
[ny, nx] = size(gx);
[Ap_RE, As_RE, Aa, Ab, As] = matrices_A_intrate_smoothing(param, setp, hx);


%% min FEVs by choice of gain gbar

tic
%Optimization Parameters
options = optimset('fmincon');
options = optimset(options, 'TolFun', 1e-9, 'display', 'none');

gbar0 = gbar; % start at CEMP's value
gbar_opt = nan(N,1);
ub = 0.2;
lb = 0.00001;
for n=1:1%N
    if mod(n,100)==0
    disp([num2str(n) ,' out of ', num2str(N)])
    end
    % Sequence of innovations
    e = [squeeze(eN(:,:,n)); zeros(1,T)]; % adding zero shocks to interest rate lag
    
    %Compute the objective function one time with some values
%     loss = obj_minFEV(gbar,gx,hx,SIG,T,burnin,e,Aa,Ab,As);
    %Declare a function handle for optimization problem
    objh = @(gbar) obj_minFEV(gbar,gx,hx,SIG,T,burnin,e,Aa,Ab,As); % solve objective for variables given params
    gbar_opt(n) = fmincon(objh, gbar0, [],[],[],[],lb,ub,[],options);
end
nanmean(gbar_opt)
toc



