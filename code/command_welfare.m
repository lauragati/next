% command_welfare.m
% compare loss in optimal Taylor rule, Taylor rule with psi_pi = 1.5 and
% optimal policy, all for the anchoring model.
% 30 July 2020

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


% fegrid_uneven = x{1};
% fegrid = fegrid_uneven;
% If you wanna use the uniform grid, then uncomment the following 3 lines:
% fegrid = linspace(femin,femax,nfe);
% x = cell(1,1);
% x{1} = fegrid;

% evaluate gradients of estimated LOM gain beforehand
k1_opt = ndim_simplex_eval(x,fegrid_fine,alph_opt);
g_fe = gradient(k1_opt);


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

% gen all the N sequences of shocks at once.
rng(0)
eN = randn(ne,T,N);
% vN = zeros(ne+1,T,N); % deal with this later




disp('%%%%%%%%%%%%%%%% Optimal PEA policy %%%%%%%%%%%%%%%%')
%% Optimal PEA policy

% Remove TR
knowTR=0
mpshock=0



%%%%%%%%%%%%%%%%%%%%%%%%%
% RE model
sig_r = param.sig_r;
sig_i = param.sig_i;
sig_u = param.sig_u;

% RE model
[fyn, fxn, fypn, fxpn] = model_NK(param);
[gx,hx]=gx_hx_alt(fyn,fxn,fypn,fxpn);
[ny, nx] = size(gx);
SIG = eye(nx).*[sig_r, sig_i, sig_u]';
eta = SIG; %just so you know

ndrop=0;

start_PEA = datetime('now');
y= zeros(ny,T,N);
for n=1:N
    disp('%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%')
    disp(['Shock sequence n=', num2str(n), ' out of N=', num2str(N)])
    disp('%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%')

    e = squeeze(eN(:,:,n));
    [~, y(:,:,n)] = sim_learnLH_pea_optimal(alph,x,param,gx,hx,SIG, PLM, gain, T,ndrop,e,g_fe,knowTR,mpshock);

end
end_PEA = datetime('now');
duration_PEA_loss = end_PEA - start_PEA

% for T=100, N=1, this takes 78 sec. 
% for T=100, N=100, it should thus take 130 min
% for T=400, N=1, it 20 minutes.
% for T=400, N=100, it should take 34 hours.

loss_PEA = loss_CB(param,y) 

