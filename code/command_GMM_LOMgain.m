% command_GMM_LOMgain
% Split up elements of materials33 to continue to estimate the approximated anchoring function
% 17 June 2020

clearvars
close all
clc

% Add all the relevant paths and grab the codename
this_code = mfilename;
[current_dir, basepath, BC_researchpath,toolpath,export_figpath,figpath,tablepath,datapath, inputsRyan_path] = add_paths;
todays_date = strrep(datestr(today), '-','_');
nowstr = strrep(strrep(strrep(datestr(now), '-','_'), ' ', '_'), ':', '_');

% Variable stuff ---
print_figs        = 0;
stop_before_plots = 0;
skip_old_plots    = 0;
output_table = print_figs;

skip = 1;
[fs, lw] = plot_configs;
redo_data_load_and_bootstrap = 0;
datestr(now)

%% Compute weighting matrix and initialize alpha
% filename ='acf_data_11_Jun_2020'; % real data

% filename ='acf_sim_data_18_Jun_2020'; % simulated data 2*15 true parameters, full Om
% filename ='acf_sim_data_19_Jun_2020'; % simulated data: 2*6 true parameters, only own Om
filename = 'acf_sim_data_21_Jun_2020'; % simulated data, 2*6 true parameters, full Om
load([filename, '.mat'])
Om = acf_outputs{1}; % this is the moments vector
Om_boot = acf_outputs{2}; % moments vectors in bootstrapped samples
ny = acf_outputs{3};
p = acf_outputs{4};
K = acf_outputs{5};
filt_data = acf_outputs{6};
lost_periods = acf_outputs{7};
T = length(filt_data);
T = T+lost_periods; % to make up for the loss due to filtering

if contains(filename,'sim')
alph_true = acf_outputs{8};
true_nk1  = acf_outputs{9};
true_nfe  = acf_outputs{10};
end


% Note: 3 moments at lag 0 are repeated. So technically we only have 42
% moments (and they could be correlated further)
% reshape(Om(1:9),3,3)

% return

% weighting matrix for GMM
% note: 2nd argument of var(X,W,DIM) signifies whether to normalize by N-1
% or N. 0 -> N-1, 1 -> N. The third argument, DIM, says along which
% dimension to take the variance.
W = diag(var(Om_boot,0,2));
W1 = W^(-1);

% return
[param, setp, param_names, param_values_str, param_titles] = parameters_next;
ne = 3;

sig_r = param.sig_r;
sig_i = param.sig_i;
sig_u = param.sig_u;

% RE model
[fyn, fxn, fypn, fxpn] = model_NK(param);
[gx,hx]=gx_hx_alt(fyn,fxn,fypn,fxpn);
[ny, nx] = size(gx);

% [Aa, Ab, As] = matrices_A_13_true_baseline(param, hx);
SIG = eye(nx).*[sig_r, sig_i, sig_u]';
eta = SIG; %just so you know

% Params for the general learning code
constant_only = 1; % learning constant only
constant_only_pi_only = 11; % learning constant only, inflation only
mean_only_PLM = -1;
slope_and_constant = 2;

dgain = 1;  % 1 = decreasing gain, 2 = endogenous gain, 3 = constant gain
again_critCEMP  = 21;
again_critCUSUM = 22;
again_critsmooth = 23;
cgain = 3;

% Model selection
%%%%%%%%%%%%%%%%%%%
PLM = constant_only_pi_only;
gain = again_critsmooth;
%%%%%%%%%%%%%%%%%%%

% display which model you're doing
[PLM_name, gain_name, gain_title] = give_names(PLM, gain);

% Specify info assumption on the Taylor rule and not to include a monpol
% shock
knowTR =1
mpshock=1
%%%%%%%%%%%%%%%%%%%

% Size of cross-section
% we're not doing a whole cross-section here
ndrop = 5 % 0-50

% gen all the N sequences of shocks at once.
rng(0)
e = randn(ne,T+ndrop); % turned monpol shocks on in smat.m to avoid stochastic singularity!

% Create grids
nk1 = 2;
nfe = 6 % 6,9,12,15
% grids for k^(-1)_{t-1} and f_{t|t-1}
k1min = 0;
k1max = 1; 
femax = 5;
femin = -femax;
k1grid = linspace(k1min,k1max,nk1); 
fegrid = linspace(femin,femax,nfe); % for alph0, fe is between (-2.6278,3.5811). N=100 AR(1)-simulations of the model yield an average fe in (-0.2946,0.2809)

% map to ndim_simplex
x = cell(2,1);
x{1} = k1grid;
x{2} = fegrid;
[xxgrid, yygrid] = meshgrid(k1grid,fegrid);

% values for k^{-1}_t for the grid
kmesh = fk_smooth_pi_only(param,yygrid,1./xxgrid); % I've checked and this gives the same as putting fe and k_{t-1} one-by-one thru fk_smooth_pi_only
k1 = 1./kmesh;

% Do an initial approx of the anchoring function to initialize the coeffs
alph0 = ndim_simplex(x,[xxgrid(:)';yygrid(:)'],k1);
rng(5)
alph0 = rand(size(alph0));
alph0 = 1.1*ones(size(alph0));

% Let's plot the approximated evolution of the gain on a finer sample
ng_fine = 100;
k1grid_fine = linspace(k1min,k1max,ng_fine);
fegrid_fine = linspace(femin,femax,ng_fine);
[xxgrid_fine, yygrid_fine] = meshgrid(k1grid_fine,fegrid_fine);

k10 = ndim_simplex_eval(x,[xxgrid_fine(:)';yygrid_fine(:)'],alph0);

xlabel = '$k^{-1}_{t-1}$'; ylabel = '$fe_{t|t-1}$'; zlabel = '$k^{-1}_{t}$';
figname = [this_code, '_initial_approx_', todays_date];
% create_pretty_3Dplot(k10,xxgrid_fine,yygrid_fine,xlabel,ylabel,zlabel,figname,print_figs)

% % % seems to behave fine for initial alphas
% [x0, y0, k0, phi0, FA0, FB0, FEt_10, diff0] = sim_learnLH_clean_approx(alph0,x,param,gx,hx,eta, PLM, gain, T,ndrop,e,knowTR,mpshock);
% % 
% % Some titles for figures
% seriesnames = {'\pi', 'x','i'};
% invgain = {'Inverse gain'};
% figname = [this_code, '_initial_obs_',PLM_name,'_', todays_date];
% create_plot_observables(y0,seriesnames, '', figname, print_figs)
% figname = [this_code, '_initial_gain_',PLM_name,'_', todays_date];
% % figname = [this_code, '_true_gain_',PLM_name,'_', todays_date];
% create_plot_observables(1./k0,invgain, '', figname, print_figs)


% return

%% GMM

% dbstop if error
% dbstop if warning

% Fmincon
%Optimization Parameters
options = optimset('lsqnonlin');
options = optimset(options, 'TolFun', 1e-9, 'display', 'iter'); 
% options.MaxFunEvals = 30000;
% options.MaxIter = 1200;
options.UseParallel = 1; % 2/3 of the time


% let's keep these bounds
ub = k1max*ones(size(alph0));
lb = k1min*ones(size(alph0));

% %Compute the objective function one time with some values
% let's weight the prior...
Wprior=0; 
[res0, Om0] = obj_GMM_LOMgain(alph0,x,xxgrid_fine,yygrid_fine,param,gx,hx,eta,e,T,ndrop,PLM,gain,p,Om,W1, alph0, Wprior);
tic
%Declare a function handle for optimization problem
objh = @(alph) obj_GMM_LOMgain(alph,x,xxgrid_fine,yygrid_fine,param,gx,hx,eta,e,T,ndrop,PLM,gain,p,Om,W1, alph0, Wprior);
[alph_opt,resnorm,residual,flag] = lsqnonlin(objh,alph0,lb,ub,options);
toc


% Let's add the final output to the finer sample
k1_opt = ndim_simplex_eval(x,[xxgrid_fine(:)';yygrid_fine(:)'],alph_opt);
disp('%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%')
disp('Is optimal k1 ever negative?')
find(k1_opt<0)

if flag==1 || flag== 2 || flag==3 % only plot if converged to a root
xlabel = '$k^{-1}_{t-1}$'; ylabel = '$fe_{t|t-1}$'; zlabel = '$k^{-1}_{t}$';
figname = [this_code, '_estimated_alphas_','ndrop_', num2str(ndrop), '_nfe_', num2str(nfe), '_', todays_date];
create_pretty_3Dplot(k1_opt,xxgrid_fine,yygrid_fine,xlabel,ylabel,zlabel,figname,print_figs)
end

% how does the model behave for estimated alpha?
[x0, y0, k0, phi0, FA0, FB0, FEt_10, diff0] = sim_learnLH_clean_approx(alph_opt,x,param,gx,hx,eta, PLM, gain, T,ndrop,rand(size(e)),knowTR,mpshock);
disp('%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%')
disp('Is implied simulated k1 ever negative?')
find(k0<0)

% Some titles for figures
seriesnames = {'\pi', 'x','i'};
invgain = {'Inverse gain'};
create_plot_observables(y0,seriesnames, 'Simulation using estimated LOM-gain approx', [this_code, '_plot1_',PLM_name,'_', todays_date], 0)
create_plot_observables(1./k0,invgain, 'Simulation using estimated LOM-gain approx', [this_code, '_plot1_',PLM_name,'_', todays_date], 0)

% Plot ACFs at start and end (Om0 and Om1 are the model-implied moments, initial and optimal)
[res1, Om1] = obj_GMM_LOMgain(alph_opt,x,xxgrid_fine,yygrid_fine,param,gx,hx,eta,e,T,ndrop,PLM,gain,p,Om,W1, alph0, Wprior);
yfig = [Om'; Om0'; Om1'];
figname= [this_code, '_ACFs_', todays_date];
create_pretty_plot_holdon(yfig,{'data', 'initial', 'optimal'},figname,print_figs)


if contains(filename,'sim')
alph_true-alph_opt
sum(abs(alph_true-alph_opt))

% plot true relationship
k1_true = ndim_simplex_eval(x,[xxgrid_fine(:)';yygrid_fine(:)'],alph_true);
% create_pretty_3Dplot(k1_true,xxgrid_fine,yygrid_fine,xlabel,ylabel,zlabel,['true_relationship', todays_date],print_figs)

% plot true, original and estimated alphas
yfig = [alph_true'; alph0'; alph_opt'];
figname= [this_code, '_alphas_', todays_date];
create_pretty_plot_holdon(yfig,{'true', 'initial', 'optimal'},figname,print_figs)

end

return

estim_outputs = {xxgrid_fine,yygrid_fine, ng_fine, k1_opt, alph_opt, x, boundname, ndrop};
filename = ['estim_LOMgain_outputs', nowstr];
save([filename, '.mat'], 'estim_outputs')
disp(['Saving as ' filename])