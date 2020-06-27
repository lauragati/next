% materials34
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
datestr(now)

do2D =0
do1D =0

nsearch = 1; % each search takes about 16 sec. So 100 should take a little under 30 min.

%% Get the data, filter them, generate data moments and bootstrap to get the weighting matrix
if skip==0
    command_acf_data
    command_acf_sim_data
    command_acf_sim_data_univariate
    %% Estimate given a dataset
    
    command_GMM_LOMgain
    
    %%  Estimate univariate anchoring function
    command_GMM_LOMgain_univariate
    
end

%% Multiple starting points, % First the 2D
if do2D==1
    
    
    filename ='acf_data_11_Jun_2020'; % real data, full Om
    %     filename = 'acf_sim_data_21_Jun_2020'; % simulated data, 2*6 true parameters, full Om
    load([filename, '.mat'])
    
    nk1=2;
    nfe=6;
    k1min = 0;
    k1max=1;
    femax = 5;
    femin = -femax;
    
    % Uniform random starting values
    b=1; a=0;
    ALPH0 = a + (b-a).*rand(nk1*nfe,nsearch);
    alph_opt = ones(nk1*nfe,nsearch);
    resnorm = zeros(1,nsearch);
    res = zeros(45,nsearch);
    Om_opt = zeros(45,nsearch);
    flag = zeros(1,nsearch);
    
    % A question: where do we gain more time: parforing the outer, or the inner
    % loop? Maybe the outer.
    tic
    parfor i=1:nsearch
        alph0 = ALPH0(:,i);
        [alph_opt(:,i), resnorm(i),res(:,i), Om_opt(:,i),flag(i)] = fun_GMM_LOMgain(acf_outputs, nk1, nfe, k1min, k1max, femin, femax, alph0);
    end
    toc
    
    % plot estimated alphas and estimation residuals
    table_flags = [(1:nsearch)', flag'];
    table_norms = [(1:nsearch)',resnorm'];
    
    searchID = num2str(1):num2str(nsearch);
    legendentries = num2cell(searchID);
    figname = [this_code,'_', filename '_2D_alphas_opt_', todays_date];
    create_pretty_plot_holdon(alph_opt(:,flag>0)',legendentries,figname,print_figs)
    
    figname = [this_code,'_', filename '_2D_residuals_', todays_date];
    create_pretty_plot_holdon(res',legendentries,figname,print_figs)
    
    % plot moments
    Om_true = acf_outputs{1};
    Om_mean = mean(Om_opt,2);
    legendentries = {'data moments', 'mean optimal moments'};
    figname = [this_code,'_', filename '_2D_moments_', todays_date];
    create_pretty_plot_holdon([Om_true';Om_mean'],legendentries,figname,print_figs)
    
    % means and medians
    alph_mean = mean(alph_opt,2);
    sorted = sort(alph_opt,2);
    if mod(nsearch,2)==0
        alph_med = (sorted(:,nsearch/2) + sorted(:,nsearch/2+1))/2;
    else
        alph_med = sorted(:,ceil(nsearch/2));
    end
    if contains(filename,'sim')
        alph_true = acf_outputs{8};
        nfe_true  = acf_outputs{9};
        figname = [this_code, '_2D_alphas_true', todays_date];
        create_pretty_plot_x(1:nk1*nfe,alph_true',figname,print_figs)
    else
        k1grid = linspace(k1min,k1max,nk1);
        fegrid = linspace(femin,femax,nfe);
        x = cell(2,1);
        x{1} = k1grid;
        x{2} = fegrid;
        [xxgrid, yygrid] = meshgrid(k1grid,fegrid);
        ng_fine = 100;
        k1grid_fine = linspace(k1min,k1max,ng_fine);
        fegrid_fine = linspace(femin,femax,ng_fine);
        [xxgrid_fine, yygrid_fine] = meshgrid(k1grid_fine,fegrid_fine);
        % optimal relationship at mean alph_opt
        k1_opt_mean = ndim_simplex_eval(x,[xxgrid_fine(:)';yygrid_fine(:)'],alph_mean);
        xlabel = '$k^{-1}_{t-1}$'; ylabel = '$fe_{t|t-1}$'; zlabel = '$k^{-1}_{t}$';
        figname = [this_code,'_', filename '_relationship_at_opt_mean_', todays_date];
        create_pretty_3Dplot(k1_opt_mean,xxgrid_fine,yygrid_fine,xlabel,ylabel,zlabel,figname,print_figs)
        % optimal relationship at median alph_opt
        k1_opt_med = ndim_simplex_eval(x,[xxgrid_fine(:)';yygrid_fine(:)'],alph_med);
        xlabel = '$k^{-1}_{t-1}$'; ylabel = '$fe_{t|t-1}$'; zlabel = '$k^{-1}_{t}$';
        figname = [this_code,'_', filename '_relationship_at_opt_median_', todays_date];
        create_pretty_3Dplot(k1_opt_med,xxgrid_fine,yygrid_fine,xlabel,ylabel,zlabel,figname,print_figs)
        
    end
    
end
%% Then the 1D
if do1D==1
    
    filename ='acf_data_11_Jun_2020'; % real data, full Om
    %      filename = 'acf_sim_univariate_data_21_Jun_2020'; % simulated data, nfe = 6. full Om
    
    load([filename, '.mat'])
    
    nfe=5;
    k1min = 0;
    k1max=1;
    femax = 3.5;
    femin = -femax;
    
    % Uniform random starting values
    b=1; a=0;
    ALPH0 = a + (b-a).*rand(nfe,nsearch);
    alph_opt = ones(nfe,nsearch);
    resnorm = zeros(1,nsearch);
    res = zeros(45,nsearch);
    Om_opt = zeros(45,nsearch);
    flag = zeros(1,nsearch);
    
   
    tic
    for i=1:nsearch
        alph0 = ALPH0(:,i);
        [alph_opt(:,i), resnorm(i),res(:,i), Om_opt(:,i),flag(i)] = fun_GMM_LOMgain_univariate(acf_outputs, nfe, k1min, k1max, femin, femax, alph0);
    end
    toc
    
    % plot estimated alphas and estimation residuals
    table_flags = [(1:nsearch)', flag'];
    table_norms = [(1:nsearch)',resnorm'];
    searchID = num2str(1):num2str(nsearch);
    legendentries = num2cell(searchID);
    figname = [this_code,'_', filename '_1D_alphas_opt_', todays_date];
    create_pretty_plot_holdon(alph_opt(:,flag>0)',legendentries,figname,print_figs)
    
    figname = [this_code,'_', filename '_1D_residuals_', todays_date];
    create_pretty_plot_holdon(res',legendentries,figname,print_figs)
    
    % plot moments
    Om_true = acf_outputs{1};
    Om_mean = mean(Om_opt,2);
    legendentries = {'data moments', 'mean optimal moments'};
    figname = [this_code,'_', filename '_1D_moments_', todays_date];
    create_pretty_plot_holdon([Om_true';Om_mean'],legendentries,figname,print_figs)
    
    %means and medians
    alph_mean = mean(alph_opt,2);
    sorted = sort(alph_opt,2);
    if mod(nsearch,2)==0
        alph_med = (sorted(:,nsearch/2) + sorted(:,nsearch/2+1))/2;
    else
        alph_med = sorted(:,ceil(nsearch/2));
    end
    
    % create grids and whatnot for plotting
    % grids for f_{t|t-1}
    fegrid = linspace(femin,femax,nfe); % for alph0, fe is between (-2.6278,3.5811).
    % map to ndim_simplex
    x = cell(1,1);
    x{1} = fegrid;
    ng_fine = 100;
    fegrid_fine = linspace(femin,femax,ng_fine);
    
    
    if contains(filename,'sim')
        alph_true = acf_outputs{8};
        nfe_true  = acf_outputs{9};
        figname = [this_code, '1D_alphas_true', todays_date];
        create_pretty_plot_x(fegrid,alph_true',figname,print_figs)
    end
end
return

%% Coax the solver to get to the right answer: for 1D case only
% --> do it on the server, as it seems that for a large enough nsearch, it can actually recover the true coeffs

if skip==0
    
    nsearch=5;
    disp(['Expected to take ', num2str(nsearch*4/60), ' minutes.'])
    filename = 'acf_sim_univariate_data_21_Jun_2020'; % simulated data, nfe = 6. full Om
    load([filename, '.mat'])
    alph_true = acf_outputs{8};
    nfe=6;
    k1min = 0;
    k1max= 1;
    femax = 5;
    femin = -femax;
    
    % Uniform random starting values
    rng('default')
    b=1; a=0;
    ALPH0 = a + (b-a).*rand(nfe,nsearch);
    alph_opt = ones(nfe,nsearch);
    resnorm = zeros(1,nsearch);
    res = zeros(45,nsearch);
    Om_opt = zeros(45,nsearch);
    flag = zeros(1,nsearch);
    
    tic
    parfor i=1:nsearch
        alph0 = ALPH0(:,i);
        [alph_opt(:,i), resnorm(i),res(:,i), Om_opt(:,i),flag(i)] = fun_GMM_LOMgain_univariate(acf_outputs, nfe, k1min, k1max, femin, femax, alph0);
    end
    toc
    
    flag
    [min_resnorm,min_idx]= min(resnorm);
    [alph_true, alph_opt(:,min_idx)]
    
    
end
%% Coax, adaptively

nsearch=5;
filename = 'acf_sim_univariate_data_21_Jun_2020'; % simulated data, nfe = 6. full Om
load([filename, '.mat'])
alph_true = acf_outputs{8};
nfe=6;
k1min = 0;
k1max= 1;
femax = 5;
femin = -femax;

% Uniform random starting values
rng('default')
b=1; a=0;
ALPH0 = a + (b-a).*rand(nfe,nsearch);
alph_opt = ones(nfe,nsearch);
resnorm = zeros(1,nsearch);
res = zeros(45,nsearch);
Om_opt = zeros(45,nsearch);
flag = zeros(1,nsearch);

best10000 =   [0.9741
    1.0000
    0.8881
    0.9325
    0.9504
    0.0000];
best10000obj = 34.5273; % with these as starting values, the first iter should reject

maxiter=4;
chain_obj = inf(1,maxiter);
chain_obj(1) = best10000obj;
chain_alph = zeros(nfe,maxiter);
chain_alph(:,1) = best10000;
iter=0;
cand_obj = 1;
disp(['Expected to take at most ', num2str(maxiter*60/60), ' minutes.'])
while cand_obj > 1e-9 && iter < maxiter
    iter =iter+1
    tic
    parfor i=1:nsearch
        alph0 = ALPH0(:,i);
        [alph_opt(:,i), resnorm(i),res(:,i), Om_opt(:,i),flag(i)] = fun_GMM_LOMgain_univariate(acf_outputs, nfe, k1min, k1max, femin, femax, alph0);
    end
    toc
    
    flag
    [cand_obj,min_idx]= min(resnorm);
    cand_alph = alph_opt(:,min_idx);
    
    % if the new proposed sol is better than the previous
    if cand_obj < chain_obj(iter)
        disp('Accepting proposal')
        chain_obj(iter+1) = cand_obj;
        chain_alph(:,iter+1) = cand_alph;
        % Perturb optimal sol nsearch times
        ALPH0 = cand_alph + 0.9*randn(size(ALPH0));
    else
        disp('Rejecting proposal')
        chain_obj(iter+1) = chain_obj(iter);
        chain_alph(:,iter+1) = chain_alph(:,iter);
        ALPH0 = chain_alph(:,iter) + 0.1*randn(size(ALPH0));
    end
    
end

chain_obj
[chain_alph, alph_true]



return
% end
%% Analyze stuff from server
filename = 'best_n10000_24_Jun_2020'; % simulated data, nfe = 6. full Om
load([filename, '.mat'])

true_vs_best = outputs{1};
nsearch  = outputs{2};
alph_opt = outputs{3};
resnorm  = outputs{4};
res      = outputs{5};
Om_opt   = outputs{6};
flag     = outputs{7}; % to check that flag seems to be always 0

alph_true = true_vs_best(:,1);
resnorm_plus = resnorm(resnorm>0);
[min_resnorm,min_idx]= min(resnorm_plus);

[alph_true, alph_opt(:,min_idx)]

