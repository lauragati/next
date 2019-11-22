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

% Time
T = 5*400 % 400
burnin = 0;

% Cross-section
N = 100 %500 size of cross-section
eN = randn(ne,T,N); % gen all the N sequences of shocks at once.

%% min FEVs by choice of gain gbar
skip_this=1;
if skip_this==0
    tic
    %Optimization Parameters
    options = optimset('fmincon');
    options = optimset(options, 'TolFun', 1e-9, 'display', 'none');
    
    gbar0 = gbar; % start at CEMP's value
    gbar_opt = nan(N,1);
    ub = 0.2;
    lb = 0.00001;
    for n=1:N
        if mod(n,100)==0
            disp([num2str(n) ,' out of ', num2str(N)])
        end
        % Sequence of innovations
        e = [squeeze(eN(:,:,n)); zeros(1,T)]; % adding zero shocks to interest rate lag
        
        %Compute the objective function one time with some values
        %     loss = obj_minFEV(gbar,gx,hx,SIG,T,burnin,e,Aa,Ab,As);
        %Declare a function handle for optimization problem
        oh = @(gbar) obj_minFEV(gbar,gx,hx,SIG,T,burnin,e,Aa,Ab,As);
        gbar_opt(n) = fmincon(oh, gbar0, [],[],[],[],lb,ub,[],options);
    end
    disp(['mean(gbar_opt)=',num2str(nanmean(gbar_opt)), '; var(gbar_opt)=', num2str(nanvar(gbar_opt))])
    toc
    
    %% plot MSE = FEV for various values of gbar
    
    ng = 100;
    gbar_vals = linspace(-0.0002, 0.0002, ng);
    mse = zeros(1,ng);
    
    e = [squeeze(eN(:,:,1)); zeros(1,T)];
    for i=1:ng
        gbar_i = gbar_vals(i);
        mse(i) = obj_minFEV(gbar_i,gx,hx,SIG,T,burnin,e,Aa,Ab,As);
    end
    
    figure
    set(gcf,'color','w'); % sets white background color
    set(gcf, 'Position', get(0, 'Screensize')); % sets the figure fullscreen
    plot(gbar_vals,mse)
end
%% try the exercise Ryan had in mind: simulate data using gbar = 0.145, and then let agents choose gbar_opt to min FEV
%Optimization Parameters
options = optimset('fmincon');
options = optimset(options, 'TolFun', 1e-9, 'display', 'none');

tic
PLM      = 1; % constant only
cgain    = 3; % constant gain
critCEMP = 1; % this doesn't matter
gbar0 = gbar;
ub = 0.2;
lb = 0.00001;
gbar_o = nan(N,1);
for n=1:N
    e = [squeeze(eN(:,:,n)); zeros(1,T)]; % adding zero shocks to interest rate lag
    %     [xsim, ysim] = sim_learn(gx,hx,SIG,T,burnin,e, Aa, Ab, As, param, PLM, cgain, critCEMP);
    
    % suppose the DGP was RE
    [xsim, ysim] = sim_model(gx,hx,SIG,T,burnin,e);
    
    oh = @(gbar) obj_minFEV_Ryan(gbar, xsim,ysim, gx, hx,T);
    gbar_o(n) = fmincon(oh, gbar0, [],[],[],[],lb,ub,[],options);
end
disp(['mean(gbar_opt)=',num2str(nanmean(gbar_o)), '; var(gbar_opt)=', num2str(nanvar(gbar_o))])

toc
ng = 100;
gbar_vals = linspace(-0.0002, 0.002, ng);
mse_Ryan = zeros(1,ng);

e = [squeeze(eN(:,:,1)); zeros(1,T)];
[xsim, ysim] = sim_learn(gx,hx,SIG,T,burnin,e, Aa, Ab, As, param, PLM, cgain, critCEMP);
for i=1:ng
    gbar_i = gbar_vals(i);
    mse_Ryan(i) = obj_minFEV_Ryan(gbar_i, xsim,ysim, gx, hx,T);
end

figure
set(gcf,'color','w'); % sets white background color
set(gcf, 'Position', get(0, 'Screensize')); % sets the figure fullscreen
plot(gbar_vals,mse_Ryan)
