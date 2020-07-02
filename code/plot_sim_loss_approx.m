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
skip_old_stuff = 1;

%% Parameters

% Load estimated LOM gain coefficents
filename = 'best_n100_29_Jun_2020'; % materials35 candidate
load([filename,'.mat'])
alph_best = output{1};
resnorm = output{2};
alph = alph_best(:,1);

% alph= [0.0674
%     0.0168
%     0
%     0.0168
%     0.0674]; % default*5

% grab the rest from materials35, part 2.5
nfe=5;
k1min = 0;
k1max= 1;
femax = 3.5;
femin = -femax;
% and from materials35, intro
fegrid = linspace(femin,femax,nfe);
x = cell(1,1);
x{1} = fegrid;

[param, setp, param_names, param_values_str, param_titles] = parameters_next;
ne = 3;

burnin = 0;
T = 400 % 400
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

%% Compute loss as a function of psi_pi and psi_x=0
if compute_loss==1
    % gen all the N sequences of shocks at once.
    rng(0)
    eN = randn(ne,T,N);
    
    % takes a little more than a min for M=30 (82 sec with parfor, 118 with for)
    tic
    disp(['Computing loss for psi_x=',num2str(param.psi_x),' and various values of psi_pi...'])
    M = 30;%30
    loss = zeros(1,M);
    loss_RE = zeros(1,M);
    psi_pi_vals = linspace(1,5,M);
    pis_x_here = 0;
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
    
    yseries=loss_RE;
    seriesnames = 'Loss RE';
    figname = [this_code, '_', 'loss','_', 'RE', '_',relevant_params, '_', date_today];
    figtitle = ['CB loss as a function of \psi_{\pi} ; ' , ' RE'];
    create_plot(xseries,yseries,seriesnames,figname,0,figtitle)
end
%% Pretty plots for draft or prezi
if print_figs==1
    figname = [this_code,'_pretty', '_', 'loss','_', gain_name, '_', PLM_name , '_', ...
        'T_', num2str(T), '_N_', num2str(N), '_burnin_', num2str(burnin),'_', ...
        relevant_params, '_date_',date_today];
    create_pretty_plot_x(xseries, loss,figname,print_figs)
    
    figname = [this_code, '_pretty', '_', 'loss_RE','_', gain_name, '_', PLM_name , '_',relevant_params, '_', date_today];
    create_pretty_plot_x(xseries, loss_RE,figname,print_figs)
    
    figname = [this_code, '_pretty', '_', 'losses','_', gain_name, '_', PLM_name , '_',relevant_params, '_', date_today];
    y = [loss_RE; loss];
    create_pretty_plot_x_holdon(xseries,y,{'RE', 'Anchoring'},figname,print_figs)
end