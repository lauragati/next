% plot_sim_loss_approx.m
% plot simulated loss
% a shortened version of grid_search.m that only does the loss-plotting part to save time
% for the case of the univariate approximated LOM gain
% this version begun in materials35
% 30 June 2020

clearvars
close all
clc

date_today = strrep(datestr(today, 'yyyy-mm-dd'), '-','_');

% % Add all the relevant paths and grab the codename
this_code = mfilename;
[current_dir, basepath, BC_researchpath,toolpath,export_figpath,figpath,tablepath,datapath] = add_paths;

% Variable stuff ---
print_figs        = 0;
stop_before_plots = 0;
skip_old_plots    = 0;
output_table = print_figs;

plot_simulated_sequence = 1;
compute_loss=1;

%% Parameters

% filename = 'estim_LOMgain_outputs_univariate_coax11_Sep_2020_15_46_40'; % materials44 candidate
filename = 'estim_LOMgain_outputs_univariate_coax15_Sep_2020_16_14_00'; % complete materials44 candidate (21 Sept draft)

% load the saved stuff
load([filename,'.mat'])
% Structure of saved file:
% estim_configs={nfe,gridspacing,femax,femin,ub,lb,Wprior,Wdiffs2,Wmid,Wmean,T,ndrop,N,eN, rngsetting};
% learn_configs = {param,PLM_name, gain_name, knowTR, mpshock};
% estim_outputs = {fegrid_fine, ng_fine, alph_opt, alph_k, ALPH0, x, estim_configs, learn_configs};
fegrid_fine = estim_outputs{1};
ng_fine     = estim_outputs{2};
alph_opt      = estim_outputs{3};
alph_k        = estim_outputs{4};
ALPH0         = estim_outputs{5};
x             = estim_outputs{6};
estim_configs = estim_outputs{7};
learn_configs = estim_outputs{8};
nfe            = estim_configs{1};
gridspacing    = estim_configs{2};
femax          = estim_configs{3};
femin          = estim_configs{4};
ub             = estim_configs{5};
lb             = estim_configs{6};
Wprior         = estim_configs{7};
Wdiffs2        = estim_configs{8};
Wmid           = estim_configs{9};
Wmean          = estim_configs{10};
T_est          = estim_configs{11};
ndrop_est      = estim_configs{12};
N_est          = estim_configs{13};
eN_est         = estim_configs{14};
rngsetting_est = estim_configs{15};
param       = learn_configs{1};
PLM_name    = learn_configs{2};
gain_name   = learn_configs{3};
knowTR_est  = learn_configs{4};
mpshock_est = learn_configs{5};

% return

fegrid = x{1};
alph = alph_opt;

% [param, setp, param_names, param_values_str, param_titles] = parameters_next;
ne = 3;

burnin = 0;
T = 100 % 400
% Size of cross-section
N = 100 %500

% Params for the general learning code
constant_only = 1; % learning constant only
constant_only_pi_only = 11; % learning constant only, inflation only (scalar learning)
slope_and_constant = 2;

dgain = 1;  % 1 = decreasing gain, 2 = endogenous gain, 3 = constant gain
again_critCEMP  = 21;
again_critCUSUM = 22;
again_smooth = 23;

cgain = 3;

% Model selection
%%%%%%%%%%%%%%%%%%%
PLM = constant_only_pi_only;
gain = again_smooth;
%%%%%%%%%%%%%%%%%%%
[PLM_name, gain_name, gain_title] = give_names(PLM, gain);

disp(['(psi_x, lam_x, lam_i)=   ', num2str([param.psi_x, param.lamx, param.lami])])

knowTR=1

% vN=zeros(ne+1,T,N);

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% %% 23 August 2020: just try the calibrated values in command_simgas.m (Materials 42)
% 
% alph = [1.0000    0.5000         0    0.5000    1.0000]'
% fegrid = [-4,-3,0,3,4]
% x{1} = fegrid;
% 
% [param, setp, param_names, param_values_str, param_titles] = parameters_next;
% 
% sig_r = 0.01;
% sig_i = 2;
% sig_u = 0.5;
% 
% eta = eye(3).*[sig_r, sig_i, sig_u]'
% 
% setp.sig_r = sig_r;
% setp.sig_i = sig_i;
% setp.sig_u = sig_u;
% setp.lamx  = 1;
% setp.lami  = 1;
% 
% % return

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% %% 27 August 2020: calibration C (Materials 43)
% 
% alph = [0.8    0.4         0    0.4    0.8]'
% fegrid = [-4,-3,0,3,4]
% x{1} = fegrid;
% 

%% Vary stuff
[param, setp, param_names, param_values_str, param_titles] = parameters_next;
% 
% setp.lamx  = 0;
% setp.lami  = 0;
% 
pis_x_here = param.psi_x;


% return

%% Compute loss as a function of psi_pi and psi_x=0
% dbstop if caught error
if compute_loss==1
    disp(datestr(now))
    % gen all the N sequences of shocks at once.
    rng(0)
    eN = randn(ne,T,N);
    
    % takes a little more than a min for M=30 (82 sec with parfor, 118 with for)
    tic
    disp(['Computing loss for psi_x=',num2str(param.psi_x),' and various values of psi_pi...'])
    M = 30;%30
    loss = zeros(1,M);
    loss_RE = zeros(1,M);
    psi_pi_vals = linspace(1,2,M);
%     pis_x_here = 0;
    parfor m=1:M
        if mod(m,10)==0
            disp(['Iteration ', num2str(m), ' out of ', num2str(M)])
        end
        psi_pi = psi_pi_vals(m);
        loss(m) = objective_CB_approx([psi_pi,pis_x_here],setp,eN,burnin,PLM,gain,alph,x, knowTR);
        loss_RE(m) = objective_CB_RE([psi_pi,pis_x_here],setp,eN,burnin);
    end
    toc
end
%% Plot loss
% create indexes for the positions of the parameters
fn = fieldnames(param);
make_index(fn)
interesting_param_names = param_names([psi_pi_idx, psi_x_idx, gbar_idx, lamx_idx, lami_idx]);
interesting_param_vals = param_values_str([psi_pi_idx, psi_x_idx, gbar_idx, lamx_idx, lami_idx]);
param_names_vals = cell(size(interesting_param_vals));
relevant_params = 'params';
for i=1:size(param_names_vals,2)
    param_names_vals{i} = [interesting_param_names{i},'_',interesting_param_vals{i}];
    relevant_params = [relevant_params, '_', param_names_vals{i}];
end
xseries=psi_pi_vals;

if print_figs==0
    yseries=loss;
    seriesnames = 'Loss';
    figname = [this_code, '_', 'loss','_', gain_name, '_', PLM_name , '_',relevant_params, '_', date_today];
    figtitle = ['CB loss as a function of \psi_{\pi} ; ' , gain_title];
    create_plot(xseries,yseries,seriesnames,figname,0,figtitle)
    
%     return
    yseries=loss_RE;
    seriesnames = 'Loss RE';
    figname = [this_code, '_', 'loss','_', 'RE', '_',relevant_params, '_', date_today];
    figtitle = ['CB loss as a function of \psi_{\pi} ; ' , ' RE'];
    create_plot(xseries,yseries,seriesnames,figname,0,figtitle)
end

% return
%% Pretty plots for draft or prezi
if print_figs==1
%     figname = [this_code,'_pretty', '_', 'loss','_', gain_name, '_', PLM_name , '_', ...
%         'T_', num2str(T), '_N_', num2str(N), '_burnin_', num2str(burnin),'_', ...
%         relevant_params, '_date_',date_today];
%     create_pretty_plot_x(xseries, loss,figname,print_figs)
%     
%     figname = [this_code, '_pretty', '_', 'loss_RE','_', gain_name, '_', PLM_name , '_',relevant_params, '_', date_today];
%     create_pretty_plot_x(xseries, loss_RE,figname,print_figs)
    
    figname = [this_code, '_pretty', '_', 'losses','_', gain_name, '_', PLM_name , '_','lamx', strrep(num2str(setp.lamx), '.','_'), '_lami', num2str(setp.lami), '_', date_today];
    loss_opt = 5*ones(size(loss));
    y = [loss; loss_RE; loss_opt];
    % 
    xlplus = [0.5,0.0002];
    ylplus = [-100,0]; % delete this one
    create_pretty_plot_x_holdon(xseries,y,{'Anchoring', 'RE', 'Optimal'},'$\psi_{\pi}$','Loss',xlplus, ylplus, figname,print_figs)
end