% command_GMM_LOMgain_univariate
% Same as command_GMM_LOMgain, except estimates a univariate anchoring
% function (gain specified in levels, not changes)
% 21 June 2020

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

save_estim_outputs =0;

skip = 1;
[fs, lw] = plot_configs;
redo_data_load_and_bootstrap = 0;
datestr(now)

%% Compute weighting matrix and initialize alpha
% filename ='acf_data_11_Jun_2020'; % real data
filename = 'acf_sim_univariate_data_21_Jun_2020'; % simulated data, nfe = 6. Note: I'm using the large moments vector.

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
    nfe_true  = acf_outputs{9};
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

% Create grid
%%%%%%%%%%%%%%%%%%%
nfe = 6 % 6,9,12,15
%%%%%%%%%%%%%%%%%%%
% grids for f_{t|t-1}
femax = 5;
femin = -femax;
fegrid = linspace(femin,femax,nfe); % for alph0, fe is between (-2.6278,3.5811).

% map to ndim_simplex
x = cell(1,1);
x{1} = fegrid;
[xxgrid] = meshgrid(fegrid);

% return
% values for k^{-1}_t for the grid
param.rho_k =0;
kmesh = fk_smooth_pi_only(param,xxgrid,rand(size(xxgrid))); % I've checked and this gives the same as putting fe and k_{t-1} one-by-one thru fk_smooth_pi_only
k1 = 1./kmesh;

% Do an initial approx of the anchoring function to initialize the coeffs
alph0 = ndim_simplex(x,xxgrid(:)',k1);
rng(8)
alph0 = rand(size(alph0));
% alph0 = 0.2*ones(size(alph0));
figname = [this_code, '_initial_alphas_', todays_date];
if skip==0
    create_pretty_plot_x(fegrid,alph0',figname,print_figs)
    pause(3)
    close
end

% Let's plot the initial approximated evolution of the gain on a finer sample
ng_fine = 100;
fegrid_fine = linspace(femin,femax,ng_fine);
k10 = ndim_simplex_eval(x,fegrid_fine,alph0);

if skip==0
    % simulation given initial alphas
    [x0, y0, k0, phi0, FA0, FB0, FEt_10, diff0] = sim_learnLH_clean_approx_univariate(alph0,x,param,gx,hx,eta, PLM, gain, T,ndrop,e,knowTR,mpshock);
    %
    % Some titles for figures
    seriesnames = {'\pi', 'x','i'};
    invgain = {'Inverse gain'};
    figname = [this_code, '_initial_obs_',PLM_name,'_', todays_date];
    create_plot_observables(y0,seriesnames, '', figname, print_figs)
    pause(2)
    close
    
    figname = [this_code, '_initial_gain_',PLM_name,'_', todays_date];
    create_plot_observables(1./k0,invgain, '', figname, print_figs)
    
    pause(3)
    close
    
end
% return

%% GMM

% dbstop if error
% dbstop if warning

% Fmincon
% %Optimization Parameters
% options = optimset('lsqnonlin');
% options = optimset(options, 'TolFun', 1e-9, 'display', 'iter');
% % options.MaxFunEvals = 30000;
% % options.MaxIter = 1200;
% % options.TolX = 1e-9;
% options.UseParallel = 1; % 2/3 of the time

%Optimization Parameters
options = optimoptions('lsqnonlin');
options = optimoptions(options, 'display', 'iter');
options.TolFun= 1e-11;
% options.OptimalityTolerance = 1e-9; % this is the guy you can access in optimoptions, not in optimset. It pertains to first order optimality measure.
% options.MaxFunEvals = 1000;
% options.MaxIter = 1200;
options.TolX = 1e-11;
options.UseParallel = 1; % 2/3 of the time


% let's keep these bounds
ub = ones(size(alph0));
lb = zeros(size(alph0));

% %Compute the objective function one time with some values
% let's weight the prior...
Wprior=0;
[res0, Om0] = obj_GMM_LOMgain_univariate(alph0,x,fegrid_fine,param,gx,hx,eta,e,T,ndrop,PLM,gain,p,Om,W1, alph0, Wprior);


tic
%Declare a function handle for optimization problem
objh = @(alph) obj_GMM_LOMgain_univariate(alph,x,fegrid_fine,param,gx,hx,eta,e,T,ndrop,PLM,gain,p,Om,W1, alph0, Wprior);
[alph_opt,resnorm,residual,flag] = lsqnonlin(objh,alph0,lb,ub,options);
toc



% Let's add the final output to the finer sample
k1_opt = ndim_simplex_eval(x,fegrid_fine(:)',alph_opt);
disp('%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%')
disp('Is optimal k1 ever negative?')
find(k1_opt<0)

if flag==1 || flag== 2 || flag==3 % only plot if converged to a root
    figname = [this_code, '_alph_opt_', todays_date];
    create_pretty_plot_x(fegrid,alph_opt',figname,print_figs)
end

flag

% how does the model behave for estimated alpha?
[x0, y0, k0, phi0, FA0, FB0, FEt_10, diff0] = sim_learnLH_clean_approx_univariate(alph_opt,x,param,gx,hx,eta, PLM, gain, T,ndrop,rand(size(e)),knowTR,mpshock);
disp('%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%')
disp('Is implied simulated k1 ever negative?')
find(k0<0)

% Some titles for figures
seriesnames = {'\pi', 'x','i'};
invgain = {'Inverse gain'};
if skip==0
    create_plot_observables(y0,seriesnames, 'Simulation using estimated LOM-gain approx', [this_code, '_plot1_',PLM_name,'_', todays_date], 0)
    create_plot_observables(1./k0,invgain, 'Simulation using estimated LOM-gain approx', [this_code, '_plot1_',PLM_name,'_', todays_date], 0)
end

% Plot ACFs at start and end (Om0 and Om1 are the model-implied moments, initial and optimal)
[res1, Om1] = obj_GMM_LOMgain_univariate(alph_opt,x,fegrid_fine,param,gx,hx,eta,e,T,ndrop,PLM,gain,p,Om,W1, alph0, Wprior);
yfig = [Om'; Om0'; Om1'];
figname= [this_code, '_ACFs_', todays_date];
create_pretty_plot_holdon(yfig,{'data', 'initial', 'optimal'},figname,print_figs)


if contains(filename,'sim')
    alph_true-alph_opt
    sum(abs(alph_true-alph_opt))
    
    % plot true, original and estimated alphas
    yfig = [alph_true'; alph0'; alph_opt'];
    figname= [this_code, '_alphas_', todays_date];
    create_pretty_plot_holdon(yfig,{'true', 'initial', 'optimal'},figname,print_figs)
end

% investigate loss function
% 1. loss(true coeffs)=0?
[res_true, Om_true] = obj_GMM_LOMgain_univariate(alph_true,x,fegrid_fine,param,gx,hx,eta,e,T,ndrop,PLM,gain,p,Om,W1, alph0, Wprior);
% res_true all are zero, as Peter said that they better be.
% Om - Om_true % these are also all zeros. So at least the code is ok.
% 2. What does the loss look like?
nrange =100;
alphi_values =linspace(lb(1),ub(1),nrange);
obj = zeros(length(alph_true),nrange);
tic
for i=1:length(alph_true)
    alph = alph_true;
    for j=1:nrange
        alph(i) = alphi_values(j);
        res = obj_GMM_LOMgain_univariate(alph,x,fegrid_fine,param,gx,hx,eta,e,T,ndrop,PLM,gain,p,Om,W1, alph0, Wprior);
        obj(i,j) = sum(res.^2);
    end
end
toc

[min_obj, min_idx] = min(obj,[],2);
[alph_true,alphi_values(min_idx)']

figure
set(gcf,'color','w'); % sets white background color
set(gcf, 'Position', get(0, 'Screensize')); % sets the figure fullscreen
for i=1:length(alph_true)
    subplot(2,3,i)
    plot(alphi_values,obj(i,:), 'linewidth', lw);
    ax = gca; % current axes
    ax.FontSize = fs;
    set(gca,'TickLabelInterpreter', 'latex');
    grid on
    grid minor
end
figname = [this_code,'_loss_for_indi_alphas_others_at_true', todays_date];
if print_figs ==1
    disp(figname)
    cd(figpath)
    export_fig(figname)
    cd(current_dir)
    close
end

% Can I find the right alphas if I start nearly at the correct alph0?
alph0 = alph_true;
alph0(6) = alphi_values(end); % <--- the answer is "depends where you start''
options.TolFun= 1e-11;
% options.OptimalityTolerance = 1e-9;
% options.MaxFunEvals = 1000;
% options.MaxIter = 1200;
options.TolX = 1e-11;
tic
[alph_opt,resnorm,residual,flag] = lsqnonlin(objh,alph0,lb,ub,options);
toc
[alph_true, alph_opt]


% save estimation outputs
if save_estim_outputs==1
    estim_outputs = {xxgrid_fine,yygrid_fine, ng_fine, k1_opt, alph_opt, x, boundname, ndrop};
    filename = ['estim_LOMgain_outputs_univariate', nowstr];
    save([filename, '.mat'], 'estim_outputs')
    disp(['Saving as ' filename])
end