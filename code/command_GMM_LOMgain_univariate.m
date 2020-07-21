% command_GMM_LOMgain_univariate
% Same as command_GMM_LOMgain, except estimates a univariate anchoring
% function (gain specified in levels, not changes)
% 21 June 2020

clearvars
close all
clc

% Add all the relevant paths and grab the codename
this_code = mfilename;
[current_dir, basepath, BC_researchpath,toolpath,export_figpath,figpath,tablepath,datapath] = add_paths;
todays_date = strrep(datestr(today), '-','_');
nowstr = strrep(strrep(strrep(datestr(now), '-','_'), ' ', '_'), ':', '_');

% Variable stuff ---
print_figs        = 0;
if contains(current_dir, 'gsfs0') % sirius server
    print_figs=1;
end
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
% % % % % % filename = 'acf_sim_univariate_data_21_Jun_2020'; % simulated data, nfe = 6. Note: I'm using the large moments vector.
% % % % % % filename = 'acf_sim_univariate_data_24_Jun_2020'; % simulated data, nfe=6, convex true function, alphas between 0 and 0.1.
% % % % % % filename = 'acf_sim_univariate_data_25_Jun_2020'; % simulated data, nfe=6, convex true function, alphas between 0 and 0.1, fe in (-3.5,3.5).
% filename = 'acf_sim_univariate_data_04_Jul_2020'; % simulated data, nfe=6, convex true function, alphas between 0 and 0.1, fe in (-3.5,3.5), new parameters, rng(0)
filename = 'acf_sim_univariate_data_06_Jul_2020'; % simulated data, nfe=5, fe=(-2,2), alph_true = (0.05; 0.025; 0; 0.025; 0.05); see Notes 6 July 2020
% filename = 'acf_sim_univariate_data_mean_21_Jul_2020'; % simulated data, nfe=5, fe=(-2,2), alph_true = (0.05; 0.025; 0; 0.025; 0.05); moments generated as average of 100 simulated datasets from true params

%%%%%%%%%%%%%%%%%%%
% Grid
nfe = 5 % 5,7,9
gridspacing = 'uniform'; % uniform or uneven
% grids for fe_{t|t-1}
femax = 2; % 3.5
femin = -2;
% upper and lower bounds for estimated coefficients
ub = ones(nfe,1); %1
lb = zeros(nfe,1); %0
% weights on additional moments
Wprior=0;%0
Wdiffs2= 100000;%10000000=10M, seems like 100K is sufficient, or even 10K
Wmid =0; %1000
Wmean=0;%100, 0
% rng(8)
% alph0 = rand(nfe,1);
% alph0 = 0.1*ones(nfe,1);
use_smart_alph0=1;% default
% alph0 =     [0.0674
%     0.0168
%          0
%     0.0168
%     0.0674]; % default*5
cross_section = 'Nestimations'; % Nestimations or Nsimulations

%Optimization Parameters
options = optimoptions('lsqnonlin');
options = optimoptions(options, 'display', 'none');
% options.TolFun= 1e-9;
% options.OptimalityTolerance = 1e-9; % this is the guy you can access in optimoptions, not in optimset. It pertains to first order optimality measure.
options.MaxFunEvals = 700;
% options.MaxIter = 1200;
% options.TolX = 1e-9;
options.UseParallel = 0; % 2/3 of the time
%%%%%%%%%%%%%%%%%%%

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
% T = 2*T

ndrop = 5 % 0-50

% gen all the N sequences of shocks at once.
rng(1) % rng('default')=rng(0)is the one that was used to generate the true data.
% Size of cross-section
N=100
eN = randn(ny,T+ndrop,N);

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

% Specify info assumption on the Taylor rule and whether to include a monpol
% shock
knowTR =1
mpshock=1
%%%%%%%%%%%%%%%%%%%
% turned monpol shocks on in smat.m to avoid stochastic singularity!

switch gridspacing
    case 'uniform'
        fegrid = linspace(femin,femax,nfe); % for alph0, fe is between (-2.6278,3.5811).
    case 'uneven'
        fegrid = uneven_grid(femin,femax,nfe)
end
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
if use_smart_alph0==1
    alph0 = ndim_simplex(x,xxgrid(:)',k1);
end


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
    [x0, y0, k0, phi0, FA0, FB0, FEt_10, diff0] = sim_learnLH_clean_approx_univariate(alph_true,x,param,gx,hx,eta, PLM, gain, T,ndrop,e,knowTR,mpshock);
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
% dbstop if error % with the catch block, you don't actually stop at the
% caught error


e0 = squeeze(eN(:,:,1));
% %Compute the objective function one time with some values
[res0, Om0] = obj_GMM_LOMgain_univariate(alph0,x,fegrid_fine,param,gx,hx,eta,e0,T,ndrop,PLM,gain,p,Om,W1,Wdiffs2,Wmid,Wmean);
disp(['Truth at e(:,:,1) has a residual of ', num2str(sum(res0.^2))])

% return
alph_opt = zeros(nfe,N);
resnorm  = zeros(1,N);
residual = zeros(length(res0),N);
flag     = zeros(1,N);
res1     = nan(size(residual));
Om1      = nan(length(Om0),N);
FE       = zeros(ny,T,N);


switch cross_section
    case 'Nestimations'
        tic
        parfor n=1:N
            e_n = squeeze(eN(:,:,n));
            %Declare a function handle for optimization problem
            objh = @(alph) obj_GMM_LOMgain_univariate(alph,x,fegrid_fine,param,gx,hx,eta,e_n,T,ndrop,PLM,gain,p,Om,W1,Wdiffs2,Wmid,Wmean);
            try
                [alph_opt(:,n),resnorm(n),residual(:,n),flag(n)] = lsqnonlin(objh,alph0,lb,ub,options);
            catch err
                disp(['History n = ', num2str(n)])
                fprintf(1,'The identifier was:\n%s',err.identifier);
                fprintf(1,'\n The error message was:\n%s',err.message);
                fprintf(1,'\n');
                alph_opt(:,n) = nan(nfe,1);
                resnorm(n) = inf;
                flag(n) = nan;
                continue % Pass control to the next iteration of FOR or WHILE loop.
            end
            
            if flag(n) >0
                % Om0 and Om1 are the model-implied moments, initial and optimal
                [res1(:,n), Om1(:,n), FE(:,:,n)] = obj_GMM_LOMgain_univariate(alph_opt(:,n),x,fegrid_fine,param,gx,hx,eta,e_n,T,ndrop,PLM,gain,p,Om,W1,0,Wmid,Wmean);
            end
        end
        toc
        
        flag
        
        alph_opt_conv = alph_opt(:,flag>0);
        resmean = nanmean(res1,2);
        resnorm_mean = sum(resmean.^2);
        
        alph_opt_mean = mean(alph_opt_conv,2)
        alph_sorted = sort(alph_opt_conv,2);
        med_idx = size(alph_sorted,2)/2;
        if isinteger(med_idx)
            alph_opt_med = (alph_sorted(:,med_idx) + alph_sorted(:,med_idx+1))/2
            alph_opt_med_unsorted = (alph_opt_conv(:,med_idx) +alph_opt_conv(:,med_idx+1))/2
            
        else
            alph_opt_med = alph_sorted(:,ceil(med_idx))
            alph_opt_med_unsorted = alph_opt_conv(:,ceil(med_idx))
            
        end
        min(resnorm_mean)
        
        plot_fe_histograms(FE);
        figname = [this_code, '_histogram_FE_','nfe_', num2str(nfe), '_loss_', num2str(floor(min(resnorm_mean))),...
            '_gridspacing_', gridspacing, '_Wdiffs2_', num2str(Wdiffs2),'_Wmid_', num2str(Wmid), '_', todays_date];
        if contains(filename,'sim')==1
            figname = [this_code, '_histogram_FE_sim_','nfe_', num2str(nfe), '_loss_', num2str(floor(min(resnorm_mean))),...
                '_gridspacing_', gridspacing, '_Wdiffs2_', num2str(Wdiffs2),'_Wmid_', num2str(Wmid), '_', todays_date];
        end
        if print_figs ==1
            disp(figname)
            cd(figpath)
            export_fig(figname)
            cd(current_dir)
            close
        end
        
        Om1mean = nanmean(Om1,2);
        
    case 'Nsimulations'
        
        %Declare a function handle for optimization problem
        objh = @(alph) obj_GMM_LOMgain_univariate_mean(alph,x,fegrid_fine,param,gx,hx,eta,eN,T,ndrop,PLM,gain,p,Om,W1,Wdiffs2,Wmid,Wmean,N);
        tic
        [alph_opt,resnorm,residual,flag] = lsqnonlin(objh,alph0,lb,ub,options);
        toc
        [res1, Om1, FE, Om_n] = obj_GMM_LOMgain_univariate_mean(alph0,x,fegrid_fine,param,gx,hx,eta,eN,T,ndrop,PLM,gain,p,Om,W1,0,Wmid,Wmean,N);
        
        % treat the single outcomes as the mean outcomes
        alph_opt_mean = alph_opt;
        Om1mean = Om1;
        resnorm_mean = resnorm;
        
end

% Let's add the final output to the finer sample
k1_opt = ndim_simplex_eval(x,fegrid_fine(:)',alph_opt_mean);
disp('%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%')
disp('Is optimal k1 ever negative?')
find(k1_opt<0)

% if flag==1 || flag== 2 || flag==3 % only plot if converged to a root
figname = [this_code, '_alph_opt_','loss_', num2str(floor(min(resnorm_mean))),'_nfe_',num2str(nfe), '_' todays_date];
create_pretty_plot_x(fegrid,alph_opt_mean',figname,print_figs)
% end


% how does the model behave for estimated alpha?
[x0, y0, k0, phi0, FA0, FB0, FEt_10, diff0] = sim_learnLH_clean_approx_univariate(alph_opt_mean,x,param,gx,hx,eta, PLM, gain, T,ndrop,rand(size(e0)),knowTR,mpshock);
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

% Covariogram
Gamj = reshape(Om,ny,ny,K+1);
Gamj0 = reshape(Om0,ny,ny,K+1);
Gamj1 = reshape(Om1mean,ny,ny,K+1);
cvgram = zeros(ny,K+1,ny);
cvgram0 = zeros(ny,K+1,ny);
cvgram1 = zeros(ny,K+1,ny);

titles = {'$\pi_t$', '$x_t$', '$i_t$'};
titles_k = {'$\pi_{t-k}$', '$x_{t-k}$', '$i_{t-k}$'};
it=0;
figure
set(gcf,'color','w'); % sets white background color
set(gcf, 'Position', get(0, 'Screensize')); % sets the figure fullscreen
for i=1:ny
    for j=1:ny
        it=it+1;
        sp(it)=subplot(ny,ny,it);
        pos_sp = get(sp(it), 'position');
        set(sp(it), 'position', [1, 1, 0.95, 0.95].*pos_sp );
        
        z = plot(0:K,zeros(1,K+1), 'k--', 'linewidth',lw); hold on
        h = plot(0:K,squeeze(Gamj(i,j,:)), 'linewidth', lw);
        h0 = plot(0:K,squeeze(Gamj0(i,j,:)), 'linewidth', lw);
        h1 = plot(0:K,squeeze(Gamj1(i,j,:)), 'linewidth', lw);
        ax = gca; % current axes
        ax.FontSize = fs*3/4;
        set(gca,'TickLabelInterpreter', 'latex');
        grid on
        grid minor
        title([titles{i}, ' vs. ', titles_k{j}],'interpreter', 'latex', 'fontsize', fs*3/4)
    end
end
% To avoid an undefined property 'NumColumns' on the server:
if contains(current_dir, 'BC_Research') % local
    lh = legend([h,h0,h1],{'Data', 'Initial','Optimal'},'interpreter', 'latex','Position',[0.45 -0.05 0.1 0.2], 'NumColumns',3, 'Box', 'off');
elseif contains(current_dir, 'gsfs0') % sirius server
    lh = legend([h,h0,h1],{'Data', 'Initial','Optimal'},'interpreter', 'latex','Position',[0.45 -0.05 0.1 0.2], 'Box', 'off');
end
% Note position: left, bottom, width, height
figname = [this_code, '_autocovariogram_','nfe_', num2str(nfe), '_loss_', num2str(floor(min(resnorm_mean))),...
    '_gridspacing_', gridspacing, '_Wdiffs2_', num2str(Wdiffs2),'_Wmid_', num2str(Wmid), '_', cross_section, '_', todays_date];
if contains(filename,'sim')==1
    figname = [this_code, '_autocovariogram_sim_','nfe_', num2str(nfe), '_loss_', num2str(floor(min(resnorm_mean))),...
        '_gridspacing_', gridspacing, '_Wdiffs2_', num2str(Wdiffs2),'_Wmid_', num2str(Wmid),'_', cross_section, '_', todays_date];
end
if print_figs ==1
    disp(figname)
    cd(figpath)
    export_fig(figname)
    cd(current_dir)
    close
end

% True alphas against estimated ones
if contains(filename,'sim')
    if length(alph_true)==length(alph_opt_mean)
        [alph_true,alph_opt_mean]
    else
        alph_true
        alph_opt_mean
    end
    
    figname= [this_code, '_alphas_','loss_', num2str(floor(min(resnorm_mean))),'_nfe_',num2str(nfe),...
        '_gridspacing_', gridspacing, '_Wdiffs2_', num2str(Wdiffs2),'_Wmid_', num2str(Wmid), '_', cross_section, '_', todays_date];
    
    figure
    set(gcf,'color','w'); % sets white background color
    set(gcf, 'Position', get(0, 'Screensize')); % sets the figure fullscreen
    plot(linspace(femin,femax,length(alph_true)),alph_true, 'linewidth',lw); hold on
    plot(fegrid, alph0, 'linewidth',lw)
    plot(fegrid, alph_opt_mean, 'linewidth',lw)
%     plot(fegrid, alph_opt_med, 'linewidth',lw)
%     plot(fegrid, alph_opt_med_unsorted, 'linewidth',lw)
    ax = gca; % current axes
    ax.FontSize = fs*3/4;
    set(gca,'TickLabelInterpreter', 'latex');
    grid on
    grid minor
    legend({'true', 'initial', 'estimated'}, 'location', 'southoutside', 'interpreter', 'latex')
    legend('boxoff')
    if print_figs ==1
        disp(figname)
        cd(figpath)
        export_fig(figname)
        cd(current_dir)
        close
    end
    
end

%%  investigate loss function
if skip==0
    
    % 1. loss(true coeffs)=0?
    [res_true, Om_true] = obj_GMM_LOMgain_univariate(alph_true,x,fegrid_fine,param,gx,hx,eta,e,T,ndrop,PLM,gain,p,Om,W1,Wdiffs2,Wmid,Wmean);
    % res_true all are zero, as Peter said that they better be.
    % Om - Om_true % these are also all zeros. So at least the code is ok.
    % 2. What does the loss look like?
    nrange =100;
    alphi_values =linspace(lb(1),ub(1),nrange);
    obj = zeros(length(alph_true),nrange);
    tic
    for i=1:length(alph_true)
        alph = alph_true;
        %     alph = alph0;
        %     alph = alph_opt;
        for j=1:nrange
            alph(i) = alphi_values(j);
            res = obj_GMM_LOMgain_univariate(alph,x,fegrid_fine,param,gx,hx,eta,e,T,ndrop,PLM,gain,p,Om,W1,Wdiffs2,Wmid,Wmean);
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
    % figname = [this_code,'_loss_for_indi_alphas_others_at_initial', todays_date];
    % figname = [this_code,'_loss_for_indi_alphas_others_at_alph_opt', todays_date];
    
    if print_figs ==1
        disp(figname)
        cd(figpath)
        export_fig(figname)
        cd(current_dir)
        close
    end
    
    return
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
    
end
%% save estimation outputs
if save_estim_outputs==1
    rngsetting=rng;
    estim_configs={nfe,gridspacing,femax,femin,ub,lb,Wprior,Wdiffs2,Wmid,Wmean,T,ndrop,N,eN, rngsetting};
    learn_configs = {param,PLM_name, gain_name, knowTR, mpshock};
    estim_outputs = {fegrid_fine, ng_fine, k1_opt, alph_opt_mean, x, estim_configs, learn_configs};
    filename = ['estim_LOMgain_outputs_univariate', nowstr];
    save([filename, '.mat'], 'estim_outputs')
    disp(['Saving as ' filename])
end




