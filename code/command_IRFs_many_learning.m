% command_IRFs_many_learning
% Do IRFs for many learning models (for the LH expectations side-project,
% but not only)
% adapted from command_sim_many_learning.m
% 6 Dec 2019
% Update 20 Jan 2020: this code has become a test for different extensions
% of the LH learning model (w/ vector learning tho).
% Update 21 Jan 2020:
% If you want to generate IRFs for the baseline NK learning model w/o bells
% and whistles, select 'true_baseline' as an extension and 'no_info_ass' as info
% assumption. For learning, you select 'default_learning'.

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
doEE    = 0;
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
p11 = param.p11;
p22 = param.p22;
ne = 3;

% Model selection and informational assumption
%%%%%%%%%%%%%%%%
info_ass= 'suboptimal_fcst'; % 'no_info_ass', 'myopic' (old), 'suboptimal_fcst' (materials12f), 'optimal_fcst' (materials12g) 'dont_know_TR' (materials12i)
extension = 'pil'; % 'Epi', 'pil', 'il', 'baseline' or 'true_baseline' (% true_baseline is the baseline where nx=3, not 4 with rho=0)
% or 'indexation' (baseline w/ indexation in NKPC), or 'Epi_CB',
% 'Markov_switchingTR_true_baseline', 'pill'
learning = 'default_learning'; %'default_learning', 'learn_hx', 'VARlearn'
%%%%%%%%%%%%%%%%

% Params for the general learning code
constant_only = 1; % learning constant only
mean_only_PLM = -1;
slope_and_constant = 2;
% lets alternate between these
PLM = constant_only;

dgain = 1;  % 1 = decreasing gain, 2 = endogenous gain, 3 = constant gain
again = 2;
cgain = 3;
gain = cgain;

T = 400 % 400
% Size of cross-section
N = 100 %500
dt_vals = 25; % time of imposing innovation
h = 10; % h-period IRFs


% Check what you mean by "baseline":
% if rho = 0, baseline really is baseline
% if rho > 0, baseline is the old intrate-smoothing (myopic info ass).
% true_baseline is the baseline where nx=3, not 4 with rho=0
if strcmp(extension,'baseline') && rho==0
    disp('Really doing baseline model, rho=0.')
elseif strcmp(extension,'baseline') && rho ~= 0
    disp('Really doing myopic interest rate smoothing model, rho > 0.')
    extension = 'old_intrate_smoothing';
end

% Number of states
if strcmp(extension,'pil') || strcmp(extension,'il') || strcmp(extension,'old_intrate_smoothing') || strcmp(extension,'baseline') ...
        || strcmp(extension,'indexation')
    nx=4;
elseif strcmp(extension,'Epi') || strcmp(extension, 'true_baseline') || strcmp(extension, 'Epi_CB') || strcmp(extension, 'Markov_switchingTR_true_baseline')
    nx=3;
elseif strcmp(extension, 'pill')
    nx=5;
else
    error('Unclear which model, check nx!')
end


if nx==3
    SIG = eye(nx).*[sig_r, sig_i, sig_u]';
elseif nx==4
    SIG = eye(nx).*[sig_r, sig_i, sig_u, 0]';
elseif nx==5
    SIG = eye(nx).*[sig_r, sig_i, sig_u, 0, 0]';
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

% %% RE model
%
% % Standard model with lag of interest rate in TR
% [fyn, fxn, fypn, fxpn] = model_NK_intrate_smoothing(param);
% [gx,hx]=gx_hx_alt(fyn,fxn,fypn,fxpn);
% [ny, nx] = size(gx);
%
% % param.psi_pi = 1; % cheat --> psi_pi < 1 makes expectations stable but
% % observables unstable
% % Original
% [~, ~, Aa, Ab, As] = matrices_A_intrate_smoothing(param, hx);
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

if strcmp(learning,'default_learning')
    %% Old info approach for comparison
    if strcmp(info_ass,'myopic')==1
        if strcmp(extension,'baseline') || strcmp(extension,'old_intrate_smoothing')
            % Standard model with lag of interest rate in TR
            [fyn, fxn, fypn, fxpn] = model_NK_intrate_smoothing(param);
            [gx,hx]=gx_hx_alt(fyn,fxn,fypn,fxpn);
            [ny, nx] = size(gx);
            [~, ~, Aa, Ab, As] = matrices_A_intrate_smoothing(param, hx);
            
        elseif strcmp(extension,'il')
            [fyn, fxn, fypn, fxpn] = model_NK_intrate_smoothing(param);
            [gx,hx]=gx_hx_alt(fyn,fxn,fypn,fxpn);
            [ny, nx] = size(gx);
            [~, ~, Aa, Ab, As] = matrices_A_intrate_smoothing(param, hx);
            
        elseif strcmp(extension,'pil')
            [fyn, fxn, fypn, fxpn] = model_NK_pilTR(param);
            [gx,hx]=gx_hx_alt(fyn,fxn,fypn,fxpn);
            [ny, nx] = size(gx);
            [Aa, Ab, As] = matrices_A_12h2(param, hx);
            
        elseif strcmp(extension,'indexation')
            [fyn, fxn, fypn, fxpn] = model_NK_indexation(param);
            [gx,hx]=gx_hx_alt(fyn,fxn,fypn,fxpn);
            [ny, nx] = size(gx);
            [Aa, Ab, As] = matrices_A_13indexation_myopic(param, hx);
        end
        
        
    elseif strcmp(info_ass,'suboptimal_fcst')
        %% New January 2020 matrices, "suboptimal forecasters" info assumption
        
        if strcmp(extension,'pil')
            % pil-model (lagged inflation in Taylor rule)
            [fyn, fxn, fypn, fxpn] = model_NK_pilTR(param);
            [gx,hx]=gx_hx_alt(fyn,fxn,fypn,fxpn);
            [ny, nx] = size(gx);
            [Aa, Ab, As] = matrices_A_12f2(param, hx);
        elseif strcmp(extension,'il')
            % il-model (interest rate smoothing)
            [fyn, fxn, fypn, fxpn] = model_NK_intrate_smoothing(param);
            [gx,hx]=gx_hx_alt(fyn,fxn,fypn,fxpn);
            [ny, nx] = size(gx);
            [Aa, Ab, As] = matrices_A_12f3(param, hx);
        elseif strcmp(extension,'indexation')
            [fyn, fxn, fypn, fxpn] = model_NK_indexation(param);
            [gx,hx]=gx_hx_alt(fyn,fxn,fypn,fxpn);
            [ny, nx] = size(gx);
            [Aa, Ab, As] = matrices_A_13indexation_subopt(param, hx);
        elseif strcmp(extension,'pill')
            [fyn, fxn, fypn, fxpn] = model_NK_pill(param);
            [gx,hx]=gx_hx_alt(fyn,fxn,fypn,fxpn);
            [ny, nx] = size(gx);
            [Aa, Ab, As] = matrices_A_13_pill_subopt_fcst(param, hx);
            
        end
        
    elseif strcmp(info_ass,'optimal_fcst')
        %% New January 2020 matrices, "optimal forecasters" info assumption
        
        if strcmp(extension,'pil')
            % pil-model (lagged inflation in Taylor rule)
            [fyn, fxn, fypn, fxpn] = model_NK_pilTR(param);
            [gx,hx]=gx_hx_alt(fyn,fxn,fypn,fxpn);
            [ny, nx] = size(gx);
            [Aa, Ab, As] = matrices_A_12g2(param, hx);
        elseif strcmp(extension,'il')
            % il-model (interest rate smoothing)
            [fyn, fxn, fypn, fxpn] = model_NK_intrate_smoothing(param);
            [gx,hx]=gx_hx_alt(fyn,fxn,fypn,fxpn);
            [ny, nx] = size(gx);
            [Aa, Ab, As] = matrices_A_12g3(param, hx);
        elseif strcmp(extension,'indexation')
            [fyn, fxn, fypn, fxpn] = model_NK_indexation(param);
            [gx,hx]=gx_hx_alt(fyn,fxn,fypn,fxpn);
            [ny, nx] = size(gx);
            [Aa, Ab, As] = matrices_A_13indexation_opt(param, hx);
            
        end
    elseif strcmp(info_ass,'no_info_ass')
        
        % Finally: the plain vanilla baseline
        if strcmp(extension,'true_baseline')
            [fyn, fxn, fypn, fxpn] = model_NK(param);
            [gx,hx]=gx_hx_alt(fyn,fxn,fypn,fxpn);
            [ny, nx] = size(gx);
            [Aa, Ab, As] = matrices_A_13_true_baseline(param, hx);
            
        elseif strcmp(extension,'Epi')
            % Epi-model (expected inflation in Taylor rule) - this one is the same
            % regardless of info assumption
            [fyn, fxn, fypn, fxpn] = model_NK_EpiTR(param);
            [gx,hx]=gx_hx_alt(fyn,fxn,fypn,fxpn);
            [ny, nx] = size(gx);
            [Aa, Ab, As, Ae] = matrices_A_12f1(param, hx); % Note: you don't actually
            %         need Ae because it's just [0,0,psi_pi] anyway, so I've let the learning code
            %         sim_learnLH_12f1.m simply use param.psi_pi for it.
        elseif strcmp(extension,'Epi_CB')
            % Epi_CB-model (CB's expected inflation in Taylor rule) - this one is the same
            % regardless of info assumption
            [fyn, fxn, fypn, fxpn] = model_NK_EpiTR(param);
            [gx,hx]=gx_hx_alt(fyn,fxn,fypn,fxpn);
            [ny, nx] = size(gx);
            [Aa, Ab, As, Ae] = matrices_A_12f1(param, hx); % Note: you don't actually
            %         need Ae because it's just [0,0,psi_pi] anyway, so I've let the learning code
            %         sim_learnLH_Epi_CB.m simply use param.psi_pi for it.
            %         The matrices of the Epi-extension are still valid
            %         here.
        elseif strcmp(extension, 'Markov_switchingTR_true_baseline')
            [fyn, fxn, fypn, fxpn] = model_NK_Markov_switchingTR(param);
            [gx,hx]=gx_hx_alt(fyn,fxn,fypn,fxpn);
            [ny, nx] = size(gx);
            ny = ny/2; % 2 regimes leads to 2*ny observables
            [Aa1, Ab1, As1, Aa2, Ab2, As2] = matrices_A_13_Markov_switching_true_baseline(param, hx);
            rng(0)
            r = generate_regime_sequence(p11,p22,T); % use this regime sequence for everything
            r1 = ones(1,T); % an only active regime
            r2 = 2*ones(1,T); % an only passive regime
            
            r = r2; % choose which regime sequence to use
        end
        
    elseif strcmp(info_ass,'dont_know_TR')
        % so far only baseline case is implemented ("true_baseline")
        [fyn, fxn, fypn, fxpn] = model_NK(param);
        [gx,hx]=gx_hx_alt(fyn,fxn,fypn,fxpn);
        [ny, nx] = size(gx);
        [Aa, Ab, As] = matrices_A_12i(param, hx);
    else
        error('Don''t know which info assumption requested.')
    end
    
    disp([learning, '! Model version: ', extension,'; ' ,info_ass])
    
elseif strcmp(learning,'learn_hx')
    
    if strcmp(info_ass,'no_info_ass')
        if strcmp(extension, 'true_baseline')
            % Learn hx, baseline
            [fyn, fxn, fypn, fxpn] = model_NK(param);
            [gx,hx]=gx_hx_alt(fyn,fxn,fypn,fxpn);
            [ny, nx] = size(gx);
            [Aa, Ab, As, Ba, Bb] = matrices_A_13_learnhx_baseline(param,hx);
        end
    end
    
    disp([learning, '! Model version: ', extension,'; ' ,info_ass])
    
elseif strcmp(learning, 'VARlearn')
    
    if strcmp(info_ass,'no_info_ass')
        if strcmp(extension, 'true_baseline')
            % Learn hx, baseline
            [fyn, fxn, fypn, fxpn] = model_NK(param);
            [gx,hx]=gx_hx_alt(fyn,fxn,fypn,fxpn);
            [ny, nx] = size(gx);
            [Aa, Ab, As, Ba, Bb] = matrices_A_13_learnhx_baseline(param,hx);
        end
    end
    
    disp([learning, '! Model version: ', extension,'; ' ,info_ass])
end
%% Simulate models

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
        
        % Differentiate between learning models
        if strcmp(learning,'default_learning')
            if strcmp(extension,'Epi')==0 && strcmp(extension,'Epi_CB')==0 && strcmp(extension,'Markov_switchingTR_true_baseline')==0
                % LH learning (learning both slope and constant of a vector)
                [~, y_LH] = sim_learnLH(gx,hx,SIG,T,burnin,e, Aa, Ab, As, param, PLM, gain);
            elseif strcmp(extension,'Epi')
                % Epi-version, needs a separate learning code
                [~, y_LH] = sim_learnLH_12f1(gx,hx,SIG,T,burnin,e, Aa,Ab, As, param, PLM, gain);
            elseif strcmp(extension,'Epi_CB')
                % Epi_CB-version, needs yet another separate learning code
                [~, y_LH] = sim_learnLH_Epi_CB(gx,hx,SIG,T,burnin,e, Aa,Ab, As, param, PLM, gain);
            elseif strcmp(extension, 'Markov_switchingTR_true_baseline')
                [~, y_LH] = sim_learnLH_Markov_switchingTR(gx,hx,SIG,T,burnin,e,Aa1, Ab1, As1, Aa2, Ab2, As2, r, param, PLM, gain);
            else
                warning('Model selection wasn''t clear and should''ve thrown an error before.')
            end
        elseif strcmp(learning,'learn_hx')
            [~, y_LH] = sim_learnLH_learnhx(gx,hx,SIG,T,burnin,e, Aa, Ab, As,Ba, Bb, param, PLM, gain);
        elseif strcmp(learning, 'VARlearn')
            [~, y_LH] = sim_learnLH_VARlearn(gx,hx,SIG,T,burnin,e, Aa, Ab, As, param, PLM, gain);
        else
                warning('Learning selection wasn''t clear and should''ve thrown an error before.')
        end
        
        
        if doEE==1
            % Euler equation learning (learning both slope and constant). Only this
            % is from materials3.m.
            [~, y_EE] = sim_learn_EE(gx,hx,fxpn,fxn,fypn,fyn,SIG,T,burnin,e,param, PLM, gain);
        end
        
        % Shocked
        
        % RE
        % make RE shock the same scale as learning:
        x0RE = (SIG*x0')';
        if strcmp(extension, 'Markov_switchingTR_true_baseline')==0
            [IR, iry, irx]=ir(gx,hx,x0RE,h);
            iry = iry';
            RE_fcsts = gx*hx*irx';
        elseif strcmp(extension, 'Markov_switchingTR_true_baseline')
            [IR, iry, irx]=ir_RE_Markov_switching(gx,hx,r,x0RE,h);
            iry = iry';
            % haven't implemented RE fcsts for the markov switching world
            % b/c they are state contingent
        end
        % Learning
        for t=1:nd
            dt = dt_vals(t);
            
            % Shocked
            % Differentiate between learning models
            if strcmp(learning,'default_learning')
                if strcmp(extension,'Epi')==0 && strcmp(extension,'Epi_CB')==0 && strcmp(extension,'Markov_switchingTR_true_baseline')==0
                    [~, ys_LH] = sim_learnLH(gx,hx,SIG,T,burnin,e, Aa, Ab, As, param, PLM, gain, dt, x0);
                elseif strcmp(extension,'Epi')==1
                    % Epi-version
                    [~, ys_LH] = sim_learnLH_12f1(gx,hx,SIG,T,burnin,e, Aa, Ab,As, param, PLM, gain, dt, x0);
                elseif strcmp(extension,'Epi_CB')
                    % Epi_CB-version, needs yet another separate learning code
                    [~, ys_LH] = sim_learnLH_Epi_CB(gx,hx,SIG,T,burnin,e, Aa,Ab, As, param, PLM, gain, dt, x0);
                elseif strcmp(extension, 'Markov_switchingTR_true_baseline')
                    [~, ys_LH] = sim_learnLH_Markov_switchingTR(gx,hx,SIG,T,burnin,e, Aa1, Ab1, As1, Aa2, Ab2, As2, r, param, PLM, gain, dt, x0);
                else
                    warning('Model selection wasn''t clear and should''ve thrown an error before.')
                end
            elseif strcmp(learning,'learn_hx')
                [~, ys_LH] = sim_learnLH_learnhx(gx,hx,SIG,T,burnin,e, Aa, Ab, As,Ba, Bb, param, PLM, gain, dt, x0);
            elseif strcmp(learning, 'VARlearn')
                [~, ys_LH] = sim_learnLH_VARlearn(gx,hx,SIG,T,burnin,e, Aa, Ab, As, param, PLM, gain, dt, x0);
            else
                warning('Learning selection wasn''t clear and should''ve thrown an error before.')
                
            end
            
            if doEE==1
                [~, ys_EE] = sim_learn_EE(gx,hx,fxpn,fxn,fypn,fyn,SIG,T,burnin,e,param, PLM, gain, dt, x0);
            end
            % Construct GIRs
            if doEE==1
                GIR_Y_EE(:,:,n,t) = ys_EE(:,dt:dt+h-1) - y_EE(:,dt:dt+h-1);
            end
            GIR_Y_LH(:,:,n,t) = ys_LH(:,dt:dt+h-1) - y_LH(:,dt:dt+h-1);
        end
        
    end
end
% warning on
% Construct RIRs by simple method: means (Option 1)
if doEE==1
    RIR_Y_EE = squeeze(nanmean(GIR_Y_EE,3));
end
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
    if doEE==1
        % 1) OBSERVABLES EE against RE
        series(1,:,:) = RIR_Y_EE(:,:,t)';
        series(2,:,:) = iry';
        figname = [this_code, '_', 'RIR_EE_' shocknames{s}, PLM_name, gain_name, '_gbar_', gbar_val];
        subplot_names = titles_obs;
        legendnames = {'Learning', 'RE'};
        figtitle = ['EE, shock imposed at t=', num2str(dt_vals(t))];
        create_subplot(series,subplot_names,figname,print_figs, figtitle, legendnames)
    end
    
    % 2) OBSERVABLES LH against RE
    series(1,:,:) = RIR_Y_LH(:,:,t)';
    series(2,:,:) = iry';
    %     figname = [this_code, '_', 'RIR_LH_' shocknames{s}, PLM_name, gain_name, '_gbar_', gbar_val];
    figname = [this_code, '_', 'RIR_LH_' shocknames{s}, '_', gain_name, '_gbar_', gbar_val, ...
        '_', learning,'_', extension, '_', info_ass, '_', PLM_name ];
    if strcmp(extension, 'Markov_switchingTR_true_baseline')
        if min(r==r1)
            regime='active';
        elseif min(r==r2)
            regime='passive';
        else regime='mixed_regime';
        end
        figname = [figname, '_', regime];
    end
    subplot_names = titles_obs;
    legendnames = {'Learning', 'RE'};
    figtitle = ['LH, shock imposed at t=', num2str(dt_vals(t))];
    create_subplot(series,subplot_names,figname,print_figs, figtitle, legendnames)
    
end
