% command_IRFs_many_learning
% Do IRFs for many learning models (for the LH expectations side-project,
% but not only)
% adapted from command_sim_many_learning.m
% 6 Dec 2019

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

%% Parameters
tic
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
% P = eye(ne).*[rho_r, rho_i, rho_u]';

%%%%%%%%%%%%%%%%
info_ass= 'optimal_fcst'; % 'myopic' (old), 'suboptimal_fcst' (materials12f), 'optimal_fcst' (materials12g)
extension = 'il'; % Epi, pil or il or baseline
%%%%%%%%%%%%%%%%

% These are the possible combos:
% 1. myopic & baseline
% 2. suboptimal & Epi
% 3. suboptimal & pil
% 4. suboptimal & il
% 5. optimal & Epi = 2.
% 6. optimal & pil
% 7. optimal & il


% Based on info ass, reset extension depending on what is implemented:
if strcmp(info_ass,'myopic')
    extension = 'baseline';
elseif  strcmp(info_ass,'suboptimal_fcst') && strcmp(extension,'baseline')
    error('This combination is not implemented.')
end
% Note: baseline is the old intrate-smoothing.


if strcmp(extension,'pil')==1 || strcmp(extension,'il')==1 || strcmp(extension,'baseline')==1
    nx=4;
elseif strcmp(extension,'Epi')==1 
    nx=3;
else
    error('Unclear which model, check nx!')
end


if nx==3
    SIG = eye(3).*[sig_r, sig_i, sig_u]';
elseif nx==4
    SIG = eye(nx).*[sig_r, sig_i, sig_u, 0]';
else
    error(['nx = ', num2str(nx),', adjust SIG!'])
end

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
gbar_val_raw = num2str(gbar);
gbar_val = replace(gbar_val_raw,'.','_');

current_param_values = [rho, rho_i, alph, kapp, psi_pi, sig, gbar]
current_param_names = ['\rho', '\rho_i', '\alpha', '\kappa', '\psi_{\pi}', '\sigma', '\bar{g}']

%% RE model

% Standard model with lag of interest rate in TR
[fyn, fxn, fypn, fxpn] = model_NK_intrate_smoothing(param);
[gx,hx]=gx_hx_alt(fyn,fxn,fypn,fxpn);
[ny, nx] = size(gx);

% param.psi_pi = 1; % cheat --> psi_pi < 1 makes expectations stable but
% observables unstable
% Original
[~, ~, Aa, Ab, As] = matrices_A_intrate_smoothing(param, hx);
% % % Alternative, general A-matrices
% [Aa2, Ab2, As2] = matrices_A_intrate_smoothing2(param, hx);
% [Aa3, Ab3, As3] = matrices_A_intrate_smoothing3(param, hx);
% [Aa4, Ab4, As4] = matrices_A_12f3(param, hx);
% return
% Aa - Aa2
% Ab - Ab2
% As - As2
% % % these aren't equal ever, and they shouldn't either because this is the wrong
% % % PQ method, not using condition (*)
% Aa - Aa3
% Ab - Ab3
% As - As3
% % % they equal - perfect - they should b/c 3 uses the MN method. (even if
% % rho !=0)
% % % But they are likely still not yet the correct intrate-smoothing model (il)
% Aa - Aa4
% Ab - Ab4
% As - As4
% % actually, these last ones, done with (presumably) correct PQ equal the
% % original ones as long as rho=0, but no longer if rho unequal to 0.
% % Actually that makes sense to me because they should embody a different
% % informational assumption than the MN method.
% return
%
% % % Model with E(pi) in TR instead of pi
% % [fyn, fxn, fypn, fxpn] = model_NK_EpiTR(param);
% % [gx,hx]=gx_hx_alt(fyn,fxn,fypn,fxpn);
% % [ny, nx] = size(gx);
% [Aa, Ab, As] = matrices_A_EpiTR(param, hx);
% [Aa2, Ab2, As2, Ae2] = matrices_A_12f1(param, hx); % with the new Epi-matrices, they are clearly not equal, so the old was indeed wrong.

%% Old info approach for comparison
if strcmp(info_ass,'myopic')==1
    % reset extension to baseline b/c that's the only one implemented
    extension = 'baseline';
    % Standard model with lag of interest rate in TR
    [fyn, fxn, fypn, fxpn] = model_NK_intrate_smoothing(param);
    [gx,hx]=gx_hx_alt(fyn,fxn,fypn,fxpn);
    [ny, nx] = size(gx);
    [~, ~, Aa, Ab, As] = matrices_A_intrate_smoothing(param, hx);
    
    
elseif strcmp(info_ass,'suboptimal_fcst')==1
    %% New January 2020 matrices, "suboptimal forecasters" info assumption
    
    if strcmp(extension,'pil')==1
        % pil-model (lagged inflation in Taylor rule)
        [fyn, fxn, fypn, fxpn] = model_NK_pilTR(param);
        [gx,hx]=gx_hx_alt(fyn,fxn,fypn,fxpn);
        [ny, nx] = size(gx);
        [Aa, Ab, As] = matrices_A_12f2(param, hx);
    elseif strcmp(extension,'il')==1
        % il-model (interest rate smoothing)
        [fyn, fxn, fypn, fxpn] = model_NK_intrate_smoothing(param);
        [gx,hx]=gx_hx_alt(fyn,fxn,fypn,fxpn);
        [ny, nx] = size(gx);
        [Aa, Ab, As] = matrices_A_12f3(param, hx);
    elseif strcmp(extension,'Epi')==1
        % Epi-model (expected inflation in Taylor rule)
        [fyn, fxn, fypn, fxpn] = model_NK_EpiTR(param);
        [gx,hx]=gx_hx_alt(fyn,fxn,fypn,fxpn);
        [ny, nx] = size(gx);
        [Aa, Ab, As, Ae] = matrices_A_12f1(param, hx); % Note: you don't actually
%         need Ae because it's just psi_pi anyway, so I've let the learning code
%         sim_learnLH_12f1.m simply use param.psi_pi for it.
    end
    
elseif strcmp(info_ass,'optimal_fcst')==1
    %% New January 2020 matrices, "optimal forecasters" info assumption
  
    if strcmp(extension,'pil')==1
        % pil-model (lagged inflation in Taylor rule)
        [fyn, fxn, fypn, fxpn] = model_NK_pilTR(param);
        [gx,hx]=gx_hx_alt(fyn,fxn,fypn,fxpn);
        [ny, nx] = size(gx);
        [Aa, Ab, As] = matrices_A_12g2(param, hx);
    elseif strcmp(extension,'il')==1
        % il-model (interest rate smoothing)
        [fyn, fxn, fypn, fxpn] = model_NK_intrate_smoothing(param);
        [gx,hx]=gx_hx_alt(fyn,fxn,fypn,fxpn);
        [ny, nx] = size(gx);
        [Aa, Ab, As] = matrices_A_12g3(param, hx);
    elseif strcmp(extension,'Epi')==1
        % Epi-model (expected inflation in Taylor rule) - this one is the same
    % regardless of info assumption
        [fyn, fxn, fypn, fxpn] = model_NK_EpiTR(param);
        [gx,hx]=gx_hx_alt(fyn,fxn,fypn,fxpn);
        [ny, nx] = size(gx);
        [Aa, Ab, As, Ae] = matrices_A_12f1(param, hx); % Note: you don't actually
%         need Ae because it's just psi_pi anyway, so I've let the learning code
%         sim_learnLH_12f1.m simply use param.psi_pi for it.
    end
else
    error('Don''t know which info assumption requested.')
end

disp(['Model version: ', extension,'; ' ,info_ass])

%% Simulate models

% Params for the general learning code
constant_only = 1; % learning constant only
mean_only_PLM = -1;
slope_and_constant = 2;
% lets alternate between these
PLM = constant_only;
if  PLM == constant_only
    PLM_name = 'constant_only'
elseif PLM == mean_only_PLM
    PLM_name = 'mean_only_PLM';
elseif PLM == slope_and_constant
    PLM_name = 'slope_and_constant'
end

dgain = 1;  % 1 = decreasing gain, 2 = endogenous gain, 3 = constant gain
again = 2;
cgain = 3;
gain = cgain;
if  gain == dgain
    gain_name = 'dgain'
elseif gain == cgain
    gain_name = 'cgain'
end


% % critCEMP=1;
% % critCUSUM=2;

T = 400 % 400
% Size of cross-section
N = 100 %500
eN = randn(ne,T,N); % gen all the N sequences of shocks at once.

dt_vals = 25; % time of imposing innovation
nd = size(dt_vals,2);
d = 1; % the innovation, delta
h = 10; % h-period IRFs
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
        else
            warning('nx = ', num2str(nx),', adjust shocks e!')
        end
        % Unshocked
        dbstop if warning
        if strcmp(extension,'Epi')==0
            % LH learning (learning both slope and constant of a vector)
            [~, y_LH] = sim_learnLH(gx,hx,SIG,T,burnin,e, Aa, Ab, As, param, PLM, gain);
        elseif strcmp(extension,'Epi')==1
            % Epi-version, needs a separate learning code
            [~, y_LH] = sim_learnLH_12f1(gx,hx,SIG,T,burnin,e, Aa,Ab, As, param, PLM, gain);
        else
            warning('Model selection wasn''t clear and should''ve thrown an error before.')
        end
        % Euler equation learning (learning both slope and constant). Only this
        % is from materials3.m.
        [~, y_EE] = sim_learn_EE(gx,hx,fxpn,fxn,fypn,fyn,SIG,T,burnin,e,param, PLM, gain);
        
        % Shocked
        % RE
        [IR, iry, irx]=ir(gx,hx,x0,h);
        iry = iry';
        RE_fcsts = gx*hx*irx';
        % Learning
        for t=1:nd
            dt = dt_vals(t);
            % Shocked
            if strcmp(extension,'Epi')==0
                
                [~, ys_LH] = sim_learnLH(gx,hx,SIG,T,burnin,e, Aa, Ab, As, param, PLM, gain, dt, x0);
            elseif strcmp(extension,'Epi')==1
                % Epi-version
                [~, ys_LH] = sim_learnLH_12f1(gx,hx,SIG,T,burnin,e, Aa, Ab,As, param, PLM, gain, dt, x0); %
            else
                warning('Model selection wasn''t clear and should''ve thrown an error before.')
            end
            [~, ys_EE] = sim_learn_EE(gx,hx,fxpn,fxn,fypn,fyn,SIG,T,burnin,e,param, PLM, gain, dt, x0);
            
            % Construct GIRs
            GIR_Y_EE(:,:,n,t) = ys_EE(:,dt:dt+h-1) - y_EE(:,dt:dt+h-1);
            GIR_Y_LH(:,:,n,t) = ys_LH(:,dt:dt+h-1) - y_LH(:,dt:dt+h-1);
        end
        
    end
end
% warning on
% Construct RIRs by simple method: means (Option 1)
RIR_Y_EE = squeeze(nanmean(GIR_Y_EE,3));
RIR_Y_LH = squeeze(mean(GIR_Y_LH,3));

disp(['(psi_x, psi_pi, rho, sig)=   ', num2str([psi_x, psi_pi, rho, sig])])
toc

if stop_before_plots==1
    return
end


%% Plots

shocknames = {'natrate', 'monpol','costpush'};
titles_obs = {'Inflation','Output gap','Int. rate'};
titles_fcsts = {'E^m_t(\pi_{t+1})', 'E^e_t(\pi_{t+1})'};
titles_LH = {'fa', 'fb'};
titles_FEs = {'FE^m_t(\pi_{t+1})', 'FE^e_t(\pi_{t+1})'};


for t=1:nd % for the two diff times of imposing the shock
    dt = dt_vals(t);
    
    clear series
    % 1) OBSERVABLES EE against RE
    series(1,:,:) = RIR_Y_EE(:,:,t)';
    series(2,:,:) = iry';
    figname = [this_code, '_', 'RIR_EE_' shocknames{s}, PLM_name, gain_name, '_gbar_', gbar_val];
    subplot_names = titles_obs;
    legendnames = {'Learning', 'RE'};
    figtitle = ['EE, shock imposed at t=', num2str(dt_vals(t))];
    create_subplot(series,subplot_names,figname,print_figs, figtitle, legendnames)
    
    % 2) OBSERVABLES LH against RE
    series(1,:,:) = RIR_Y_LH(:,:,t)';
    series(2,:,:) = iry';
    figname = [this_code, '_', 'RIR_LH_' shocknames{s}, PLM_name, gain_name, '_gbar_', gbar_val];
    subplot_names = titles_obs;
    legendnames = {'Learning', 'RE'};
    figtitle = ['LH, shock imposed at t=', num2str(dt_vals(t))];
    create_subplot(series,subplot_names,figname,print_figs, figtitle, legendnames)
    
end
