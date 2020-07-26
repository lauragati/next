% command_acf_sim_data_univariate_mean.m
% A version of command_acf_sim_data_univariate that simulates the true
% model N times and computes moments as a mean of that.
% - moments of the mean data, or
% - mean moments of the cross-sectional data? % <-- this feels more natural to me, so this is what I'll do
% Takes about 70 seconds altogether.
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

%% 1.) Simulate data and filter it
% Simulate data from the model given true alphas


% Some thought as to how to generate true alphas and true grids!
% Create grids
nfe = 5 % 6,9,12,15
% grids for f_{t|t-1}
femax = 2;
femin = -femax;
fegrid = linspace(femin,femax,nfe);

rng(0)
% alph_true = [0.1;0.05;0.001;0.00001;0.02;0.09];
alph_true = [0.05;0.025;0;0.025;0.05];
% alph_true = (0.005*fegrid.^2)';

% map to ndim_simplex
x = cell(1,1);
x{1} = fegrid;

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
mpshock=1  % turned monpol shocks on to avoid stochastic singularity!
%%%%%%%%%%%%%%%%%%%

% Size of cross-section
% we're not doing a whole cross-section here
ndrop = 5 % 0-50

Torig=156;
N = 100
% gen all the N sequences of shocks at once.
rng(0)
eN = randn(ny,Torig+ndrop,N);
nobs=ny+1;
vN = randn(nobs,Torig+ndrop,N);

lost_periods = 24; % 2*12
filt_data = nan(nobs,Torig-lost_periods-1,N);
Om_n = nan(nobs*nobs*5,N);
p=4; % need to fix p b/c otherwise my Bs and res's are of different sizes
B_n = nan(nobs*p+1,nobs,N); % (no. vars * no. lags + constant) by no. of vars
res_n = nan(Torig-lost_periods-p-1,nobs,N);
tic
% Takes about 4 sec for N=100
parfor n=1:N
    e = squeeze(eN(:,:,n));
    v = squeeze(vN(:,:,n));
    [x0, y0, k0, phi0, FA0, FB0, FEt_10, diff0] = sim_learnLH_clean_approx_univariate(alph_true,x,param,gx,hx,eta, PLM, gain, Torig+ndrop,ndrop,e,v,knowTR,mpshock);
    
    % return
    
    %     ng_fine = 100;
    %     fegrid_fine = linspace(femin,femax,ng_fine);
    %
    %     k10 = ndim_simplex_eval(x,fegrid_fine(:)',alph_true);
    
    
    y_data = [y0(:,1:end-1); squeeze(phi0(1,1,1:end-1))'];
    [~,T] = size(y_data)
    
    % Filter the data
    % HP filter
    g_data = nan(size(y_data));
    c_data = nan(size(y_data));
    for i=1:nobs
        [g_data(i,:),c_data(i,:)] = HPfilter(y_data(i,:)');
    end
    
    % Hamilton filter
    h=8;
    v_data = nan(nobs,T-4-h+1);
    for i=1:nobs
        [v_data(i,:)] = Hamiltonfilter(y_data(i,:)');
    end
    lost_periods_H = h+3;
    
    
    % BK filter
    K=12;
    ystar_data = nan(nobs,T-2*K);
    for i=1:nobs
        ystar_data(i,:) = BKfilter(y_data(i,:)');
    end
    lost_periods_BK = 2*K;
    
    % choose your favorite filtered data
    filt_data(:,:,n) = ystar_data;
    % lost_periods = lost_periods_BK; % manually at end
    
    
    % compute moments, Om, as the autocovariances of the data for lags
    % 0,1,...,K
    K=4;
    % Take the initial data, estimate a VAR
    %     max_lags = 16;
    %     [AIC,BIC,HQ] = aic_bic_hq(filt_data(:,:,n)',max_lags);
    %     p =min([AIC,BIC,HQ]);
    % A is the impact matrix, identified via Cholesky, B is the beta_ols, res are
    % the residuals, sigma is the estimated VC matrix.
    [A,B,res,sigma] = sr_var(filt_data(:,:,n)', p);
    
    % Rewrite the VAR(p) as VAR(1) (see Hamilton, p. 273, Mac)
    c = B(1,:); % coefficients of the constant
    PHI = B(2:end,:)';
    F = [PHI; [eye(nobs*(p-1)), zeros(nobs*(p-1),nobs)]];
    Q = [[sigma, zeros(nobs,nobs*(p-1))]; zeros(nobs*(p-1),nobs*p)];
    % check sizes
    np = nobs*p;
    sizeF = size(F) == [np,np];
    sizeQ = size(Q) == [np,np];
    
    vecSig = (eye(np^2)-kron(F,F))\vec(Q);
    % VC matrix of data y
    Gamj = zeros(nobs,nobs,K+1);
    Gamj_own = zeros(nobs,K+1);
    
    Sig = reshape(vecSig,np,np);
    for j=0:K
        % jth Autocov of data y, still as a VAR(1)
        Sigj = F^j * Sig;
        % jth Autocov of data y, finally as a VAR(p)
        Gamj(:,:,j+1) = Sigj(1:nobs,1:nobs);
        % This only takes the covariances wrt the same var: (e.g Cov(x_t,x_{t-1}))
        Gamj_own(:,j+1) = diag(Sigj(1:nobs,1:nobs));
    end
    % moments vector
    Om = vec(Gamj);
    % Om = vec(Gamj_own);
    
    
    % save for outside parfor
    B_n(:,:,n) = B;
    res_n(:,:,n) = res;
    Om_n(:,n) = Om;
    
end
toc

% Mean moments for output
Om = nanmean(Om_n,2);

% Take silly means for bootstrap
B = mean(B_n,3);
res = mean(res_n,3);

% return

%% 2.) Bootstrap the data, and get variance of moments (autocovariances from 0 to lag K) from the bootstrapped sample
% Takes 60 sec b/c parallel
% This section is inspired by main_file_SVAR_just_IT_controllingNEWS.m in my work with Marco

% Resample the residuals and use beta_ols from VAR to create nboot new samples.
nboot =10000;
q=16; % blocksize
nburn = 1000;
which_correction ='blocks';
disp(datestr(now))
disp('Creating the bootstrapped sample: takes about 60 sec')
tic
[dataset_boot] = data_boot(B,nburn,res,nboot,which_correction,q);
toc
% Autocov matrix from bootstrapped sample for lags 0,...,K
K = 4;
Gamj = zeros(nobs,nobs,K+1,nboot);
Om_boot = zeros(length(Om),nboot);
tic
disp(datestr(now))
disp('Creating the bootstrapped autocovariances: takes about 30 sec b/c parallel')
parfor i=1:nboot
    Gamj_booti = zeros(nobs,nobs,K+1);
    Gamj_own_booti = zeros(nobs,K+1);
    %         max_lags = 16;
    %         [AIC,BIC,HQ] = aic_bic_hq(squeeze(dataset_boot(:,:,i)),max_lags);
    %         p = min([AIC,BIC,HQ]); % lag selection (p) is the lag
    [A,B,res,sigma] = sr_var(squeeze(dataset_boot(:,:,i)), p);
    
    % Rewrite the VAR(p) as VAR(1) (see Hamilton, p. 273, Mac)
    c = B(1,:); % coefficients of the constant
    PHI = B(2:end,:)';
    F = [PHI; [eye(nobs*(p-1)), zeros(nobs*(p-1),nobs)]];
    v = [res'; zeros(nobs*(p-1),size(res',2))];
    Q = [[sigma, zeros(nobs,nobs*(p-1))]; zeros(nobs*(p-1),nobs*p)];
    % check sizes
    np = nobs*p;
    sizeF = size(F) == [np,np];
    sizeQ = size(Q) == [np,np];
    
    vecSig = (eye(np^2)-kron(F,F))\vec(Q);
    % VC matrix of data y
    Sig = reshape(vecSig,np,np);
    for j=0:K
        % jth Autocov of data y, still as a VAR(1)
        Sigj = F^j * Sig;
        % jth Autocov of data y, finally as a VAR(p)
        Gamj_booti(:,:,j+1) = Sigj(1:nobs,1:nobs);
        Gamj_own_booti(:,j+1) = diag(Sigj(1:nobs,1:nobs));
    end
    % gather the ACF of each bootstrapped sample
    Gamj(:,:,:,i) = Gamj_booti;
    Om_boot(:,i) = vec(Gamj_booti);
    %     Om_boot(:,i) = vec(Gamj_own_booti);
end
toc
filename = ['acf_sim_univariate_data_mean_', todays_date];
acf_outputs = {Om, Om_boot,nobs,p,K, filt_data, lost_periods, alph_true, nfe};
save([filename,'.mat'],'acf_outputs')
disp(['Saving as ' filename])