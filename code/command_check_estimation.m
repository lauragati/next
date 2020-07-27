% command_check_estimation.m
% a file to check the estimation procedure in
% comman_GMM_LOMgain_univariate.m: is it able to identify params of a
% simpler model?

% Take a stupid state-space model with a piecewise linear aprox to the
% observation equation

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
output_table = print_figs;

skip = 1;
[fs, lw] = plot_configs;
datestr(now)

%% Generate the data

rng(0)
rho=0.7;
T=100;
ndrop=5;
e = randn(1,T+ndrop);
% alph_true = [0.05;0.025;0;0.025;0.05];
alph_true = [0.025;0.05];

nknots = 2;
statemax = 2;
statemin = -statemax;
stategrid = linspace(statemin,statemax,nknots);
% map to ndim_simplex
x = cell(1,1);
x{1} = stategrid;

% states and jumps
s = zeros(T,1);
y = zeros(T,1);

s(1) = e(1);
y(1) = 0;
for t=2:T+ndrop
    s(t)=rho*s(t-1) +e(t);
    y(t)= ndim_simplex_eval(x, s(t),alph_true);
end

y = y(ndrop+1:end);

% Map into command_acf_sim_data_univariate.m
y_data = y';
[nobs,T] = size(y_data)

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
filt_data = c_data;
% lost_periods = lost_periods_BK;


% compute moments, Om, as the autocovariances of the data for lags
% 0,1,...,K
K=4;
% Take the initial data, estimate a VAR
max_lags = 16;
[AIC,BIC,HQ] = aic_bic_hq(filt_data',max_lags);
p =min([AIC,BIC,HQ]);
% A is the impact matrix, identified via Cholesky, B is the beta_ols, res are
% the residuals, sigma is the estimated VC matrix.
[A,B,res,sigma] = sr_var(filt_data', p);

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
plot(Om)

% return

% 2.) Bootstrap the data, and get variance of moments (autocovariances from 0 to lag K) from the bootstrapped sample
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
disp('Creating the bootstrapped autocovariances: takes about 80 sec b/c parallel')
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

%% Estimation
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

% dbstop if error
% dbstop if warning


W = diag(var(Om_boot,0,2));
W1 = W^(-1);

% Size of cross-section
N=100
ub = ones(nknots,1); %1
lb = zeros(nknots,1); %0
cross_section = 'Nsimulations'; % Nestimations or Nsimulations


rng(1) % rng(1)  vs. rng('default')=rng(0) is the one that was used to generate the true data.
eN = randn(1,T+ndrop,N);

alph0 = alph_true+0.01*rand(size(alph_true));

disp(['Beginning estimation...', datestr(now)])

switch cross_section
    case 'Nsimulations'
        % Compute the objective function one time with some values
        [res0, Om0, Om_n0] = obj_GMM_check_mean(alph0,x,rho,eN,T,ndrop,p,Om,W1);
        %Declare a function handle for optimization problem
        objh = @(alph) obj_GMM_check_mean(alph,x,rho,eN,T,ndrop,p,Om,W1);
        tic
        [alph_opt,resnorm,residual,flag] = lsqnonlin(objh,alph0,lb,ub,options);
        toc
        [res1, Om1, Om_n1] = obj_GMM_check_mean(alph_opt,x,rho,eN,T,ndrop,p,Om,W1);
        alph_opt_mean = alph_opt;
    case 'Nestimations'
        alph_opt = zeros(nknots,N);
        resnorm  = zeros(1,N);
        flag     = zeros(1,N);
        
        % Compute the objective function one time with some values
        [res0, Om0] = obj_GMM_check(alph0,x,rho,eN,T,ndrop,p,Om,W1);
        
        Om1      = nan(length(Om0),N);
        residual = zeros(length(res0),N);
        res1     = nan(size(residual));
        
        tic
        parfor n=1:N
            e=squeeze(eN(:,:,n));
            %Declare a function handle for optimization problem
            objh = @(alph) obj_GMM_check(alph,x,rho,e,T,ndrop,p,Om,W1);
            [alph_opt(:,n),resnorm(n),residual(:,n),flag(n)] = lsqnonlin(objh,alph0,lb,ub,options);
            
            [res1(:,n), Om1(:,n)] = obj_GMM_check(alph_opt(:,n),x,rho,e,T,ndrop,p,Om,W1);
            
        end
        toc
        alph_opt_conv = alph_opt(:,flag>0);
        resmean = nanmean(res1,2);
        resnorm_mean = sum(resmean.^2);
        
        alph_opt_mean = mean(alph_opt_conv,2);
end

flag

[alph_true, alph_opt_mean]

figure
plot(alph_true); hold on
plot(alph_opt_mean)