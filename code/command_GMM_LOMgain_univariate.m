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
investigate_loss=0;
[fs, lw] = plot_configs;
datestr(now)

%% Compute weighting matrix and initialize alpha
% filename ='acf_data_11_Jun_2020'; % real data
% filename ='acf_data_21_Jul_2020'; % real data with SPF expectation in it
% % % % % % filename = 'acf_sim_univariate_data_21_Jun_2020'; % simulated data, nfe = 6. Note: I'm using the large moments vector.
% % % % % % filename = 'acf_sim_univariate_data_24_Jun_2020'; % simulated data, nfe=6, convex true function, alphas between 0 and 0.1.
% % % % % % filename = 'acf_sim_univariate_data_25_Jun_2020'; % simulated data, nfe=6, convex true function, alphas between 0 and 0.1, fe in (-3.5,3.5).
% % % % % % filename = 'acf_sim_univariate_data_04_Jul_2020'; % simulated data, nfe=6, convex true function, alphas between 0 and 0.1, fe in (-3.5,3.5), new parameters, rng(0)
% filename = 'acf_sim_univariate_data_06_Jul_2020'; % simulated data, nfe=5, fe=(-2,2), alph_true = (0.05; 0.025; 0; 0.025; 0.05); see Notes 6 July 2020
% % % % % % filename = 'acf_sim_univariate_data_mean_21_Jul_2020'; % simulated data, nfe=5, fe=(-2,2), alph_true = (0.05; 0.025; 0; 0.025; 0.05); moments generated as average of 100 simulated datasets from true params
% filename = 'acf_sim_univariate_data_22_Jul_2020'; % simulated data with expectation in it, nfe=5, fe=(-2,2), alph_true = (0.05; 0.025; 0; 0.025; 0.05). W/ measurement error
% % % % % % filename = 'acf_sim_univariate_data_mean_26_Jul_2020'; % simulated data with expectation in it, nfe=5, fe=(-2,2), alph_true = (0.05; 0.025; 0; 0.025; 0.05); moments generated as average of 100 simulated datasets from true params
% filename = 'acf_sim_univariate_data_04_Aug_2020'; % simulated data, nfe=5, fe=(-2,2), alph_true = (0.05; 0.025; 0; 0.025; 0.05); Expectations, yes, measurement error, no!
% % % % % % filename = 'acf_sim_univariate_data_09_Aug_2020'; % simulated data, nfe=5, fe=(-0.5,0.5), alph_true = (0.05; 0.025; 0; 0.025; 0.05); Expectations, yes, measurement error, no, RIDGE.
filename = 'acf_sim_univariate_data_10_Aug_2020'; % simulated data, nfe=5, fe=(-2,2), alph_true = 4*(0.05; 0.025; 0; 0.025; 0.05); referring to point (R,b) in my notes.
% % % % % % filename = 'acf_sim_univariate_data_13_Aug_2020'; % simulated data, nfe=5, fe=(-2,2), alph_true = (0.05; 0.025; 0; 0.025; 0.05); sig_u = 2. (RIDGE)
% filename = 'acf_sim_univariate_data_18_Aug_2020'; % simulated data, nfe=5, fe=(-2,2), alph_true = (0.05; 0.025; 0; 0.025; 0.05); sig_u = 2, lam =0.001 (corrected (now standardized regressors) RIDGE)

%%%%%%%%%%%%%%%%%%%
% Grid
nfe = 5 % 5,7,9
gridspacing = 'manual'; % uniform or uneven, or manual
% grids for fe_{t|t-1}
femax = 2; % 3.5
femin = -femax;
% upper and lower bounds for estimated coefficients
ub = ones(nfe,1); %1
lb = zeros(nfe,1); %0
% weights on additional moments
Wprior=0;%0
Wdiffs2= 100000;%100000, seems like 100K is sufficient, or even 10K
Wmid =1000; %1000
Wmean=0;%100, 0
% rng(8)
% alph0 = rand(nfe,1);
% alph0 = 0.1*ones(nfe,1);
use_smart_alph0=1;% default
cross_section = 'Nsimulations'; % Nestimations or Nsimulations (default)
est_shocks = 0;
if est_shocks==1
    cross_section = 'Nsimulations';
end
scaleW =0; %0
use_expectations_data=0; %1
sig_v = 0; %0 vs 1 variance of measurement error: set to zero to shut measurement error off (default)

%Optimization Parameters
options = optimoptions('lsqnonlin');
options = optimoptions(options, 'display', 'iter');
% options.TolFun= 1e-9; % objective function tolerance, default 1.0000e-06
% options.OptimalityTolerance = 1e-9; % this is the guy you can access in
% optimoptions, not in optimset. It pertains to first order optimality
% measure. Default 1.0000e-06
options.MaxFunEvals = 700; % 700
% options.MaxIter = 33;
% options.TolX = 1e-9; % step tolerance: default 1.0000e-06
options.UseParallel = 0; % 2/3 of the time
h_sig = 100000;
h_alph= 100000;
% options.FiniteDifferenceStepSize = sqrt(eps)*[h_sig;h_sig;h_sig;h_alph;h_alph;h_alph;h_alph;h_alph;]; % default is sqrt(eps)
%%%%%%%%%%%%%%%%%%%

load([filename, '.mat'])
Om = acf_outputs{1}; % this is the moments vector
Om_boot = acf_outputs{2}; % moments vectors in bootstrapped samples
nobs = acf_outputs{3};
p = acf_outputs{4};
K = acf_outputs{5};
filt_data = acf_outputs{6};
lost_periods = acf_outputs{7};
T = length(filt_data);
T = T+lost_periods; % to make up for the loss due to filtering

ndrop = 5 % 0-50

% return

% Size of cross-section
N=100

if use_expectations_data==0 % take out moments pertaining to expectations
    Ommatrix = reshape(Om,nobs,nobs,K+1);
    % do the same for the bootstrapped variances
    nboot=size(Om_boot,2);
    Om_boot_matrix = reshape(Om_boot, nobs,nobs,K+1, nboot);
    nobs = 3;
    Ommatrix_net = Ommatrix(1:nobs,1:nobs,:);
    Om_boot_matrix_net = Om_boot_matrix(1:nobs,1:nobs,:,:);
    Om = vec(Ommatrix_net);
    [omb_m, omb_n, ~, ~] = size(Om_boot_matrix_net);
    Om_boot = reshape(Om_boot_matrix_net,omb_m*omb_n*(K+1), nboot);
    
end

% return

rng(1) % rng(1)  vs. rng('default')=rng(0) is the one that was used to generate the true data.


% gen all the N sequences of shocks at once.
eN = randn(3,T+ndrop,N);
vN = sig_v*randn(nobs,T+ndrop,N); % measurement error

% return

if contains(filename,'sim')
    alph_true = acf_outputs{8};
    nfe_true  = acf_outputs{9};
end

% Note: 3 moments at lag 0 are repeated. So technically we only have
% n_moment-3 moments (and they could be correlated further)

% return

% weighting matrix for GMM
% note: 2nd argument of var(X,W,DIM) signifies whether to normalize by N-1
% or N. 0 -> N-1, 1 -> N. The third argument, DIM, says along which
% dimension to take the variance.
sigboot = diag(var(Om_boot,0,2));
W = sigboot;

% % % Check where the not-rescaling or rescaling causes problems
% X=W;
% scaler = floor(log10(min(diag(W))));
% % a=10^(abs(scaler));
% a=10;
% Y = a*X;
% 1./Y == (1/a) ./X % aha! Not always true!
% inv(Y) == 1/a *inv(X) % but it's b/c 1/a is too close to zero.
%
% 1/a *inv(X) == inv(a*X) % not always true either! So this means that inverting the unrescaled W indeed caused troubles!
% 1/a ./X == 1./(a*X)
% inv(Y) == inv(a*X) % this is always true: so once you scale X, inverting works the way it should. But not before.
%
% inv(X)./inv(Y)
% inv(X)/inv(Y)
%
% % '/' is much more accurate than inv() or ^(-1)!
% X = [1,2; 3,4]
% inv(X) % correct
% 1./X % incorrect wtf!
% inv(X) == 1./X % aren't equal!!! Of course not. Idiot.

% If W really small, scale it up by the exponent of the smallest value
% get exponent of 10 in the smallest diagonal element
% return
if scaleW==1
    scaler = floor(log10(min(diag(sigboot))));
    if scaler < 0
        W_alt = W.* 10^(abs(scaler)); % just to check elementwise, but it gives the same thing
        W = W* 10^(abs(scaler));
    end
end
W1 = W^(-1);
% W1 = W1*10^(abs(scaler)); % scale it back up - do we get what we would w/o scaling? YES.
% W1 - inv(sigboot) % equal on the order of e-12

% % Just checking that it's true for W1 too, and it is
% W1 == 1/a *inv(X) % not always tru
% W1 == inv(a*X) % this is always true
% W1=eye(size(W1));

% Per Ryan meeting 12 August 2020:
W1 = sqrt(W1); % elementwise sqrt.

% Also per Ryan meeting 12 August 2020:
maxSig = max(max(diag(sigboot)));
minSig = min(min(diag(sigboot)));
maxminratio = maxSig/minSig; % it's < 10e+6 if no expectations are used, on the order of 9e+7 if expectations are used


% return

[param, setp, param_names, param_values_str, param_titles] = parameters_next;

sig_r = param.sig_r;
sig_i = param.sig_i;
sig_u = param.sig_u;

% RE model
[fyn, fxn, fypn, fxpn] = model_NK(param);
[gx,hx]=gx_hx_alt(fyn,fxn,fypn,fxpn);
[~, nx] = size(gx);

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
    case 'manual'
        %         fegrid = [-2,-1,1,2]
        %         fegrid = [-2,-1.5,0,1.5,2]
        fegrid = [-4,-3,0,3,4]
        %         fegrid = [-4,-3, 3,4]
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


% Finer sample
ng_fine = 100;
broaden = 2
fegrid_fine = linspace(femin-broaden,femax+broaden,ng_fine);
k10 = ndim_simplex_eval(x,fegrid_fine,alph0);

% % extrapolation works beautifully
% ndim_simplex_eval(x,[-3:1:3],alph0)


% % Can recover scaling factor?
% res0 = obj_GMM_LOMgain_univariate_mean(alph_true,x,fegrid_fine,param,gx,hx,eta,eN,vN,T,ndrop,PLM,gain,p,Om,W1,Wdiffs2,Wmid,Wmean,use_expectations_data,N);
% resnorm0 = sum(res0.^2);
% % return
% scaler = floor(log10(min(diag(sigboot))));
% sclf = 10^(abs(scaler));
% res1 = obj_GMM_LOMgain_univariate_mean(alph_true,x,fegrid_fine,param,gx,hx,eta,eN,vN,T,ndrop,PLM,gain,p,Om,sqrt(inv(W*sclf)),Wdiffs2,Wmid,Wmean,use_expectations_data,N);
% resnorm1 = sum(res1.^2);

% resnorm0 - resnorm1*sclf
% resnorm0 - resnorm1*sclf^2 % yes you do


% return

%% GMM

% dbstop if error % with the catch block, you don't actually stop at the
% caught error

% return

switch cross_section
    case 'Nestimations'
        alph_opt = zeros(nfe,N);
        resnorm  = zeros(1,N);
        flag     = zeros(1,N);
        FE       = zeros(nx,T,N);
        
        
        e0 = squeeze(eN(:,:,1));
        v0 = squeeze(vN(:,:,1));
        % Compute the objective function one time with some values
        [res0, Om0] = obj_GMM_LOMgain_univariate(alph0,x,fegrid_fine,param,gx,hx,eta,e0,v0,T,ndrop,PLM,gain,p,Om,W1,Wdiffs2,Wmid,Wmean,use_expectations_data);
        disp(['Truth at e(:,:,1) has a residual of ', num2str(sum(res0.^2))])
        
        Om1      = nan(length(Om0),N);
        residual = zeros(length(res0),N);
        res1     = nan(size(residual));
        
        
        
        tic
        parfor n=1:N
            e_n = squeeze(eN(:,:,n));
            v_n = squeeze(vN(:,:,n));
            %Declare a function handle for optimization problem
            objh = @(alph) obj_GMM_LOMgain_univariate(alph,x,fegrid_fine,param,gx,hx,eta,e_n,v_n,T,ndrop,PLM,gain,p,Om,W1,Wdiffs2,Wmid,Wmean,use_expectations_data);
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
                [res1(:,n), Om1(:,n), FE(:,:,n)] = obj_GMM_LOMgain_univariate(alph_opt(:,n),x,fegrid_fine,param,gx,hx,eta,e_n,v_n,T,ndrop,PLM,gain,p,Om,W1,Wdiffs2,Wmid,Wmean,use_expectations_data);
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
        else
            alph_opt_med = alph_sorted(:,ceil(med_idx))
            
        end
        resnorm_mean
        
        if skip==0
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
        end
        
        Om1mean = nanmean(Om1,2);
        
        
        
    case 'Nsimulations'
        if est_shocks==0
            % Testing whether moments converge to something as N--> inf, and it
            % seems yes, at N=1000.
            %         N=100; % N=10000: errors on the order of 0.005, plots almost on top of each other, max(diff)=0.0083. That seems sufficient. Is N=1000 sufficient? It seems so: plots on top, max(diff=0.0109)
            %         rng(1)
            %         eN1 = randn(3,T+ndrop,N);
            %         vN1 = randn(nobs,T+ndrop,N);
            %
            %         rng(2)
            %         eN2 = randn(3,T+ndrop,N);
            %         vN2 = randn(nobs,T+ndrop,N);
            %
            %         tic
            %         [res1, Om1, FE1, Om_n1] = obj_GMM_LOMgain_univariate_mean(alph0,x,fegrid_fine,param,gx,hx,eta,eN1,vN1,T,ndrop,PLM,gain,p,Om,W1,Wdiffs2,Wmid,Wmean,N);
            %         [res2, Om2, FE2, Om_n2] = obj_GMM_LOMgain_univariate_mean(alph0,x,fegrid_fine,param,gx,hx,eta,eN2,vN2,T,ndrop,PLM,gain,p,Om,W1,Wdiffs2,Wmid,Wmean,N);
            %         toc
            %
            %         diffOm = Om1-Om2;
            %
            %         figure
            %         plot(Om1)
            %         hold on
            %         plot(Om2)
            %         max(diffOm)
            %         return
            
            dbstop if warning
            % Compute the objective function one time with some values
            [res0, Om0, FE0, Om_n0] = obj_GMM_LOMgain_univariate_mean(alph0,x,fegrid_fine,param,gx,hx,eta,eN,vN,T,ndrop,PLM,gain,p,Om,W1,Wdiffs2,Wmid,Wmean,use_expectations_data,N);
            resnorm0 = sum(res0.^2)
            
            [res0, Om0, FE0, Om_n0] = obj_GMM_LOMgain_univariate_mean(alph_true,x,fegrid_fine,param,gx,hx,eta,eN,vN,T,ndrop,PLM,gain,p,Om,W1,Wdiffs2,Wmid,Wmean,use_expectations_data,N);
            resnorm0 = sum(res0.^2)
%             return
            %Declare a function handle for optimization problem
            objh = @(alph) obj_GMM_LOMgain_univariate_mean(alph,x,fegrid_fine,param,gx,hx,eta,eN,vN,T,ndrop,PLM,gain,p,Om,W1,Wdiffs2,Wmid,Wmean,use_expectations_data,N);
            tic
            [alph_opt,resnorm,residual,flag] = lsqnonlin(objh,alph0,lb,ub,options);
            toc
            [res1, Om1, FE, Om_n] = obj_GMM_LOMgain_univariate_mean(alph_opt,x,fegrid_fine,param,gx,hx,eta,eN,vN,T,ndrop,PLM,gain,p,Om,W1,Wdiffs2,Wmid,Wmean,use_expectations_data,N);
            
            flag
            nancount = sum(sum(isnan(Om_n)));
            nanpercent = nancount/numel(Om_n)
            % treat the single outcomes as the mean outcomes
            alph_opt_mean = alph_opt;
            Om1mean = Om1;
            resnorm_mean = resnorm
            
        elseif est_shocks==1
            disp('Estimating shock volatilities too')
            sig0 = 2*ones(nx,1);
            ub_sig = 8*ones(nx,1);
            lb_sig = 0*ones(nx,1);
            
            % call thet the full estimated parameter vector
            thet0 = [sig0; alph0]; 
            ub = [ub_sig;ub];
            lb = [lb_sig;lb];
            
            [res0, Om0, FE0, Om_n0] = obj_GMM_LOMgain_univariate_mean_shocks(thet0,x,fegrid_fine,param,gx,hx,eta,eN,vN,T,ndrop,PLM,gain,p,Om,W1,Wdiffs2,Wmid,Wmean,use_expectations_data,N);
            resnorm0 = sum(res0.^2)
            
            %Declare a function handle for optimization problem
            objh = @(thet) obj_GMM_LOMgain_univariate_mean_shocks(thet,x,fegrid_fine,param,gx,hx,eta,eN,vN,T,ndrop,PLM,gain,p,Om,W1,Wdiffs2,Wmid,Wmean,use_expectations_data,N);
            tic
            [thet_opt,resnorm,residual,flag] = lsqnonlin(objh,thet0,lb,ub,options);
            toc
            [res1, Om1, FE, Om_n] = obj_GMM_LOMgain_univariate_mean_shocks(thet_opt,x,fegrid_fine,param,gx,hx,eta,eN,vN,T,ndrop,PLM,gain,p,Om,W1,Wdiffs2,Wmid,Wmean,use_expectations_data,N);
            
            flag
            nancount = sum(sum(isnan(Om_n)));
            nanpercent = nancount/numel(Om_n)
            % treat the single outcomes as the mean outcomes
            sig_opt = thet_opt(1:3)'
            alph_opt_mean = thet_opt(4:end);
            Om1mean = Om1;
            resnorm_mean = resnorm
            
        end
        
end

%% Plots
% the appendix of each FIGNAME:
figspecs = ['N_', num2str(N),'_nfe_', num2str(nfe), ...
    '_gridspacing_', gridspacing, '_Wdiffs2_', num2str(Wdiffs2),'_Wmid_', num2str(Wmid), '_', cross_section, '_', this_code, '_', nowstr];
%'_use_expectations_', num2str(use_expectations_data), '_use_meas_error_', num2str(sig_v), ...



% Let's add the final output to the finer sample
k1_opt = ndim_simplex_eval(x,fegrid_fine(:)',alph_opt_mean);
disp('%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%')
disp('Is optimal k1 ever negative?')
find(k1_opt<0)


if skip==0
    % if flag==1 || flag== 2 || flag==3 % only plot if converged to a root
    figname = ['alph_opt_',figspecs];
    create_pretty_plot_x(fegrid,alph_opt_mean',figname,print_figs)
    % end
end

% how does the model behave for estimated alpha?
[~, y1, k11, phi1, ~, ~, FEt_11, diff1] = sim_learnLH_clean_approx_univariate(alph_opt_mean,x,param,gx,hx,eta, PLM, gain, T+ndrop,ndrop,...
    squeeze(eN(:,:,1)),sig_v*squeeze(vN(:,:,1)),knowTR,mpshock);
disp('%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%')
fe = FEt_11(1,:);

disp('Is implied simulated k1 ever negative?')
find(k11<0)

% Some titles for figures
seriesnames = {'\pi', 'x','i'};
invgain = {'Inverse gain'};
if skip==0
    create_plot_observables(y1,seriesnames, 'Simulation using estimated LOM-gain approx', ['sim_obs__alph_opt_randshocks', figspecs], 0)
    create_plot_observables(1./k11,invgain, 'Simulation using estimated LOM-gain approx', ['sim_obs__alph_opt_randshocks', figspecs], 0)
end

% Covariogram
Gamj = reshape(Om,nobs,nobs,K+1);
Gamj0 = reshape(Om0,nobs,nobs,K+1);
Gamj1 = reshape(Om1mean,nobs,nobs,K+1);
cvgram = zeros(nobs,K+1,nobs);
cvgram0 = zeros(nobs,K+1,nobs);
cvgram1 = zeros(nobs,K+1,nobs);

titles = {'$\pi_t$', '$x_t$', '$i_t$', '$E_t\pi_{t+1}$'};
titles_k = {'$\pi_{t-k}$', '$x_{t-k}$', '$i_{t-k}$', '$E_{t-k}\pi_{t-k+1}$'};
it=0;
figure
set(gcf,'color','w'); % sets white background color
set(gcf, 'Position', get(0, 'Screensize')); % sets the figure fullscreen
for i=1:nobs
    for j=1:nobs
        it=it+1;
        sp(it)=subplot(nobs,nobs,it);
        pos_sp = get(sp(it), 'position');
        set(sp(it), 'position', [1, 1, 0.95, 0.95].*pos_sp );
        
        z = plot(0:K,zeros(1,K+1), 'k--', 'linewidth',lw); hold on
        h = plot(0:K,squeeze(Gamj(i,j,:)), 'linewidth', lw);
        h0 = plot(0:K,squeeze(Gamj0(i,j,:)), 'linewidth', lw);
        h1 = plot(0:K,squeeze(Gamj1(i,j,:)), 'linewidth', lw);
        ax = gca; % current axes
        ax.FontSize = fs*3/4;
        ax.YRuler.Exponent = 0; % turns off scientific notation
        set(gca,'TickLabelInterpreter', 'latex');
        grid on
        grid minor
        title([titles{i}, ' vs. ', titles_k{j}],'interpreter', 'latex', 'fontsize', fs*2/4)
    end
end
% To avoid an undefined property 'NumColumns' on the server:
if contains(current_dir, 'BC_Research') % local
    lh = legend([h,h0,h1],{'Data', 'Initial','Optimal'},'interpreter', 'latex','Position',[0.45 -0.05 0.1 0.2], 'NumColumns',3, 'Box', 'off');
elseif contains(current_dir, 'gsfs0') % sirius server
    lh = legend([h,h0,h1],{'Data', 'Initial','Optimal'},'interpreter', 'latex','Position',[0.45 -0.05 0.1 0.2], 'Box', 'off');
end


% Note position: left, bottom, width, height
figname = ['autocovariogram_' figspecs];
if contains(filename,'sim')==1
    figname = ['autocovariogram_sim_', figspecs];
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
        disp('alph_true against mean estimate')
        [alph_true,alph_opt_mean]
    else
        disp('alph_true against mean estimate')
        alph_true
        alph_opt_mean
    end
    
    figname= ['alphas_',figspecs];
    
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

%%  investigate loss function for fixing some alphas at true values
% 28 July 2020
% takes about 90 sec.
if investigate_loss==1
    datestr(now)
    
    figspecs = ['N_', num2str(N),'_nfe_', num2str(nfe), ...
        '_gridspacing_', gridspacing, '_Wdiffs2_', num2str(Wdiffs2),'_Wmid_', num2str(Wmid), '_', cross_section, '_', this_code, '_', nowstr];
    
    % 1. loss(true coeffs)=0?
    [res_at_true, Om_at_true, FE_at_true, Om_n_at_true, expl_sim_counter] = obj_GMM_LOMgain_univariate_mean(alph_true,x,fegrid_fine,param,gx,hx,eta,eN,vN,T,ndrop,PLM,gain,p,Om,W1,Wdiffs2,Wmid,Wmean,use_expectations_data,N);
    resnorm_at_true = sum(res_at_true.^2)
    % res_true all aren't zero, as Peter said that they better be. In fact,
    % increasing in N!
    expl_percent = expl_sim_counter/N *100
    
    % 2. What does the loss look like?
    nrange =10;
    incr=0.01;%0.02, 0.001
    alphi_values =linspace(0.001,0.5,nrange);
    %     alphi_values = nan(length(alph_true),nrange);
    
    obj = nan(length(alph_true),nrange);
    tic
    for i=1:length(alph_true)
        alph = alph_true;
        %         % try to center the range tightly around the true value
        %         if i==3
        %             alphi_values(i,:) = linspace(0, alph_true(i)+incr,nrange)';
        %         else
        %             alphi_values(i,:) = linspace(alph_true(i)-incr, alph_true(i)+incr,nrange)';
        %         end
        for j=1:nrange
            alph(i) = alphi_values(j);
            try
                res = obj_GMM_LOMgain_univariate_mean(alph,x,fegrid_fine,param,gx,hx,eta,eN,vN,T,ndrop,PLM,gain,p,Om,W1,Wdiffs2,Wmid,Wmean,use_expectations_data,N);
            catch
            end
            obj(i,j) = sum(res.^2);
        end
    end
    toc
    
    figure
    set(gcf,'color','w'); % sets white background color
    set(gcf, 'Position', get(0, 'Screensize')); % sets the figure fullscreen
    for i=1:length(alph_true)
        subplot(2,3,i)
        plot(alphi_values,obj(i,:), 'linewidth', lw);
        ax = gca; % current axes
        ax.FontSize = fs;
        set(gca,'TickLabelInterpreter', 'latex');
        ax.XAxis.Exponent = 0;
        %         ax.YRuler.Exponent = 0; % turns off scientific notation
        grid on
        grid minor
    end
    sgt = sgtitle(['$\alpha^{true} = $', num2str(alph_true')]);
    sgt.FontSize =fs;
    sgt.Interpreter = 'latex';
    incr_str = strrep(num2str(incr), '.','_');
    figname = ['loss_increment_',incr_str,'_', figspecs];
    % change of name, otherwise it was too long for export_fig
    
    if print_figs ==1
        disp(figname)
        cd(figpath)
        export_fig(figname)
        cd(current_dir)
        close
    end
    
    %% Now look at individual alphas - meaning have a personalized range for each element of alpha
    datestr(now)
    nrange =10;
    alphi_values = nan(length(alph_true),nrange);
    obj_indi = nan(length(alph_true),nrange);
    
    
    tic
    for i=1:5 %1:length(alph_true)
        alph = alph_true;
        % try to center the range tightly around the true value
        if i==1 || i==5
            alphi_values(i,:) = linspace(0.045,0.055, nrange)';%linspace(0.04,0.06, nrange)'
        elseif i==2 || i==4
            alphi_values(i,:) = linspace(0.01,0.025, nrange)';%
        elseif i==3
            alphi_values(i,:) = linspace(0, 0.001,nrange)'; %
        end
        for j=1:nrange
            alph(i) = alphi_values(i,j);
            try
                res = obj_GMM_LOMgain_univariate_mean(alph,x,fegrid_fine,param,gx,hx,eta,eN,vN,T,ndrop,PLM,gain,p,Om,W1,Wdiffs2,Wmid,Wmean,use_expectations_data,N);
            catch
            end
            obj_indi(i,j) = sum(res.^2);
        end
    end
    toc
    
    figure
    set(gcf,'color','w'); % sets white background color
    set(gcf, 'Position', get(0, 'Screensize')); % sets the figure fullscreen
    for i=1:length(alph_true)
        subplot(2,3,i)
        plot(alphi_values(i,:),obj_indi(i,:), 'linewidth', lw);
        ax = gca; % current axes
        ax.FontSize = fs;
        set(gca,'TickLabelInterpreter', 'latex');
        %         ax.XAxis.Exponent = 0;
        %         ax.YRuler.Exponent = 0; % turns off scientific notation
        grid on
        grid minor
    end
    sgt = sgtitle(['$\alpha^{true} = $', num2str(alph_true')]);
    sgt.FontSize =fs;
    sgt.Interpreter = 'latex';
    
    figspecs = ['N_', num2str(N),'_nfe_', num2str(nfe), ...
        '_gridspacing_', gridspacing, '_Wdiffs2_', num2str(Wdiffs2),'_Wmid_', num2str(Wmid), '_', cross_section, '_', this_code, '_', nowstr];
    
    
    figname = ['loss_indi_nrange', num2str(nrange), '_', figspecs];
    % change of name, otherwise it was too long for export_fig
    
    if print_figs ==1
        disp(figname)
        cd(figpath)
        export_fig(figname)
        cd(current_dir)
        close
    end
    
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




