% grid_search.m
% Search the paramspace for (psi_pi, psi_x) to minimize CB's loss
% Begun in materials14.m
% 26 Jan 2020

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
do_fmincon =0;
skip_old_stuff = 1;

%% Parameters
[param, setp] = parameters_next;
ne = 3;

burnin = 0;
T = 400 % 400
% Size of cross-section
N = 100 %500

% Params for the general learning code
constant_only = 1; % learning constant only
mean_only_PLM = -1;
slope_and_constant = 2;

dgain = 1;  % 1 = decreasing gain, 2 = endogenous gain, 3 = constant gain
again_critCEMP  = 21;
again_critCUSUM = 22;
cgain = 3;

% Model selection
%%%%%%%%%%%%%%%%%%%
PLM = constant_only;
gain = again_critCUSUM;
%%%%%%%%%%%%%%%%%%%
[PLM_name, gain_name, gain_title] = give_names(PLM, gain);


%% Fmincon
% takes about 9 min
tic
% gen all the N sequences of shocks at once.
rng(0)
eN = randn(ne,T,N);

%Optimization Parameters
options = optimset('fmincon');
options = optimset(options, 'TolFun', 1e-9, 'display', 'iter');

% varp0 = [1.19,0];
% ub = [1.2,1];
% lb = [1.01,-.01];
varp0 = 1.19;
ub = 1.5;
lb = 1.01;
%Compute the objective function one time with some values
loss = objective_CB(varp0,setp,eN,burnin,PLM,gain);

%Declare a function handle for optimization problem
objh = @(varp) objective_CB(varp,setp,eN,burnin,PLM,gain);
[par_opt] = fmincon(objh, varp0, [],[],[],[],lb,ub,[],options);
% fmincon(FUN,X0,A,B,Aeq,Beq,LB,UB,NONLCON,OPTIONS)
toc

%% Compute loss as a function of psi_pi and psi_x=0
if do_fmincon==1
    % takes a little more than a min
    tic
    disp('Computing loss for psi_x=0 and various values of psi_pi...')
    M = 30;
    loss = zeros(1,M);
    psi_pi_vals = linspace(1,1.1,M);
    pis_x_here = 0;
    for m=1:M
        if mod(m,10)==0
            disp(['Iteration ', num2str(m), ' out of ', num2str(M)])
        end
        psi_pi = psi_pi_vals(m);
        loss(m) = objective_CB([psi_pi,pis_x_here],setp,eN,burnin,PLM,gain);
    end
    toc
end
%% Plot loss
% create indexes for the positions of the parameters
fn = fieldnames(param);
make_index(fn)
interesting_param_names = param_names([psi_pi_idx, psi_x_idx, gbar_idx, thetbar_idx, thettilde_idx, kap_idx, alph_cb_idx]);
interesting_param_vals = param_values_str([psi_pi_idx, psi_x_idx, gbar_idx, thetbar_idx, thettilde_idx, kap_idx, alph_cb_idx]);
param_names_vals = cell(size(interesting_param_vals));
relevant_params = 'params';
for i=1:size(param_names_vals,2)
    param_names_vals{i} = [interesting_param_names{i},'_',interesting_param_vals{i}];
    relevant_params = [relevant_params, '_', param_names_vals{i}];
end

yseries=loss;
xseries=psi_pi_vals;
seriesnames = 'Loss';
figname = [this_code, '_', 'loss','_', gain_name, '_', PLM_name , '_',relevant_params, '_', date_today];
figtitle = ['CB loss as a function of \psi_{\pi} ; ' , gain_title];
create_plot(xseries,yseries,seriesnames,figname,print_figs,figtitle)
