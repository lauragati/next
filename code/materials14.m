% materials 14
% Goals:
% 1.) maybe a last attempt at getting rid of the overshooting
% 2.) a correct way of dealing with endogenous states
% 3.) a correct way of doing the projection facility
% 4.) a correction for PLM that leads to new fa & fb
% 23 Jan 2020
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

if skip_old_stuff==0
    %% Do IRFs for extensions - same code as before
    command_IRFs_many_learning
end

%% Check the b -> gx version
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
p11 = param.p11;
p22 = param.p22;
ne = 3;

% Params for the general learning code
constant_only = 1; % learning constant only
slope_and_constant = 2;
% lets alternate between these
PLM = constant_only;

dgain = 1;
cgain = 3;
gain = cgain;

T = 400 % 400
% Size of cross-section
N = 100 %500
dt_vals = 25; % time of imposing innovation
h = 10; % h-period IRFs


% % Baseline model
% [fyn, fxn, fypn, fxpn] = model_NK(param);
% [gx,hx]=gx_hx_alt(fyn,fxn,fypn,fxpn);
% [ny, nx] = size(gx);
% [Aa, Ab, As] = matrices_A_13_true_baseline(param, hx);

% pil-model, suboptimal forecasters assumption
[fyn, fxn, fypn, fxpn] = model_NK_pilTR(param);
[gx,hx]=gx_hx_alt(fyn,fxn,fypn,fxpn);
[ny, nx] = size(gx);
[Aa, Ab, As] = matrices_A_12f2(param, hx);

if nx==3
    SIG = eye(nx).*[sig_r, sig_i, sig_u]';
elseif nx==4
    SIG = eye(nx).*[sig_r, sig_i, sig_u, 0]';
elseif nx==5
    SIG = eye(nx).*[sig_r, sig_i, sig_u, 0, 0]';
else
    error(['nx = ', num2str(nx),', adjust SIG!'])
end

% Give names
if  PLM == constant_only
    PLM_name = 'constant_only'
elseif PLM == mean_only_PLM
    PLM_name = 'mean_only_PLM';
elseif PLM == slope_and_constant
    PLM_name = 'slope_and_constant'
end

if  gain == dgain
    gain_name = 'dgain'
elseif gain == cgain
    gain_name = 'cgain'
end

% % critCEMP=1;
% % critCUSUM=2;

% gen all the N sequences of shocks at once.
eN = randn(ne,T,N);

% Preallocate
nd = size(dt_vals,2);
d = 1; % the innovation, delta
GIR_Y_EE = zeros(ny,h,N,nd);
GIR_Y_LH = zeros(ny,h,N,nd);

% warning off
for s=2  %2->zoom in on monetary policy shock
    x0 = zeros(1,nx);
    x0(s) = d;
    
    for n=1:N
        % Sequence of innovations
        if nx==3
            e = squeeze(eN(:,:,n)); % <-- EPI
        elseif nx==4
            e = [squeeze(eN(:,:,n)); zeros(1,T)]; % adding zero shocks to lagged jumps become states
        elseif nx==5
            e = [squeeze(eN(:,:,n)); zeros(1,T); zeros(1,T)]; % more padding zeros
        else
            warning('nx = ', num2str(nx),', adjust shocks e!')
        end
        
        % Unshocked
        dbstop if warning
        [xsim, y_LH] = sim_learnLH_materials14(gx,hx,SIG,T,burnin,e, Aa, Ab, As, param, PLM, gain);
        
        % Shocked
        % make RE shock the same scale as learning:
        x0RE = (SIG*x0')';
        [IR, iry, irx]=ir(gx,hx,x0RE,h);
        iry = iry';
        RE_fcsts = gx*hx*irx';
        
        % Learning
        for t=1:nd
            dt = dt_vals(t);
            % Shocked
            [~, ys_LH,evening_fcst, morning_fcst, FA, FB, FEt_1, shock, diff] = sim_learnLH_materials14(gx,hx,SIG,T,burnin,e, Aa, Ab, As, param, PLM, gain, dt, x0);
            
            % Construct GIRs
            GIR_Y_LH(:,:,n,t) = ys_LH(:,dt:dt+h-1) - y_LH(:,dt:dt+h-1);
        end
        
    end
end
% warning on
% Construct RIRs by simple method: means (Option 1)
RIR_Y_LH = squeeze(mean(GIR_Y_LH,3));

disp(['(psi_x, psi_pi, rho, sig)=   ', num2str([psi_x, psi_pi, rho, sig])])


if stop_before_plots==1
    return
end


%% Plots
date_today = strrep(datestr(today, 'yyyy-mm-dd'), '-','_');

shocknames = {'natrate', 'monpol','costpush'};
titles_obs = {'Inflation','Output gap','Int. rate'};
titles_fcsts = {'E^m_t(\pi_{t+1})', 'E^e_t(\pi_{t+1})'};
titles_LH = {'fa', 'fb'};
titles_FEs = {'FE^m_t(\pi_{t+1})', 'FE^e_t(\pi_{t+1})'};


for t=1:nd % for the two diff times of imposing the shock
    dt = dt_vals(t);
    
    clear series
    
    % 2) OBSERVABLES LH against RE
    series(1,:,:) = RIR_Y_LH(:,:,t)';
    series(2,:,:) = iry';
    %     figname = [this_code, '_', 'RIR_LH_' shocknames{s}, PLM_name, gain_name, '_gbar_', gbar_val];
    figname = [this_code, '_', 'RIR_LH_' shocknames{s}, '_', gain_name, ...
        '_', PLM_name , '_','baseline_', date_today];
    subplot_names = titles_obs;
    legendnames = {'Learning', 'RE'};
    figtitle = ['LH, shock imposed at t=', num2str(dt_vals(t))];
    create_subplot(series,subplot_names,figname,print_figs, figtitle, legendnames)
    
end