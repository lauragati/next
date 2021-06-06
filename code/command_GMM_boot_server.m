% command_GMM_boot_server.m
% Is the main file of estimating alpha on 100 bootstrapped datasets on the
% server.
% Relies somewhat on command_GMM_LOMgain_univariate.m
% 6 June 2021

clearvars
close all
clc

% Add all the relevant paths and grab the codename
this_code = mfilename;
[current_dir, basepath, BC_researchpath,toolpath,export_figpath,figpath,tablepath,datapath] = add_paths;
todays_date = strrep(datestr(today), '-','_');
nowstr = strrep(strrep(strrep(datestr(now), '-','_'), ' ', '_'), ':', '_');

% Variable stuff ---
if contains(current_dir, 'gsfs0') % sirius server
    print_figs=1;
end


datestr(now)

%% Section 1 is identical to command_GMM_LOMgain_univariate.m, except I delete comments to economize on space

%% Compute weighting matrix and initialize alpha
filename = 'acf_data_23_Aug_2020'; % real data with SPF expectation in it but qoq annualized inflation rates and expectations

% Grid
nfe = 5 % 5,7,9
gridspacing = 'manual'; % uniform or uneven, or manual
ub = ones(nfe,1); %1
lb = zeros(nfe,1); %0
% weights on additional moments
Wprior=0;%0
Wdiffs2= 100000;%100000, seems like 100K is sufficient, or even 10K
Wmid =1000; %1000
Wmean=0;%100, 0
alph0 = [0.8,0.4,0,0.4,0.8]'

use_smart_alph0=0;% default
cross_section = 'Nsimulations'; % Nestimations or Nsimulations (default)
est_shocks = 0;
if est_shocks==1
    cross_section = 'Nsimulations';
end
scaleW =0; %0
ryans_rescale = 0; % Per Ryan meeting 12 August 2020.
use_expectations_data=1; %1
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
options.UseParallel = false; % 2/3 of the time
h_sig = 100000;
h_alph= 100;%100000
% options.FiniteDifferenceStepSize = sqrt(eps)*[h_alph;h_alph;h_alph;h_alph;h_alph;]; % default is sqrt(eps); sqrt(eps)*[h_sig;h_sig;h_sig;h_alph;h_alph;h_alph;h_alph;h_alph;]
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
sigboot = diag(var(Om_boot,0,2));
W = sigboot;

if scaleW==1
    scaler = floor(log10(min(diag(sigboot))));
    if scaler < 0
        W_alt = W.* 10^(abs(scaler)); % just to check elementwise, but it gives the same thing
        W = W* 10^(abs(scaler));
    end
end
W1 = W^(-1);

% Per Ryan meeting 12 August 2020:
W1 = sqrt(W1); % elementwise sqrt.

% Also per Ryan meeting 12 August 2020:
diag_sigboot = diag(sigboot);
maxSig = max(max(diag_sigboot));
minSig = min(min(diag_sigboot));
maxminratio = maxSig/minSig; % Should be and is < 10e+6 if no expectations are used, on the order of 9e+7 or 1e+8 if expectations are used
%  is < 10e+6 if expectations are annualized too in the true data!

if ryans_rescale==1
    % Ryan's rescaling strategy
    threshold_ratio = 10e6;
    sigdiffs = maxSig ./ diag_sigboot;
    % how many too small elements do we have
    n_problems = sum(sigdiffs > threshold_ratio);
    % find the n problematic elements
    where = find(sigdiffs > threshold_ratio, n_problems, 'first');
    % fix the problematic elements of sigma to be 1/10000 the largest element
    % of sigma.
    diag_sigboot(where) = maxSig/10000;
    W = diag(diag_sigboot);
    W1 = inv(W);
end

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
        fegrid = [-4,-3,0,3,4]
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

femin = min(fegrid);
femax = max(fegrid);
% Finer sample
ng_fine = 100;
broaden = 2
fegrid_fine = linspace(femin-broaden,femax+broaden,ng_fine);
k10 = ndim_simplex_eval(x,fegrid_fine,alph0);


%% Section 2 - this is different - use first NB=100 elements of Om_boot to estimate alpha_hat.
% Just use Nsimulations
% No estimation of shock volatilities

NB = 1; %100
alpha_hat_boot = zeros(nfe,NB);
%% Here it goes: for each bootstrap moments vector

% CONT HERE WITH PARFOR
for bb = 1:NB
    Om = Om_boot(:,bb); % NEW
    
    % Compute the objective function one time with some values
    [res0, Om0, FE0, Om_n0] = obj_GMM_LOMgain_univariate_mean(alph0,x,fegrid_fine,param,gx,hx,eta,eN,vN,T,ndrop,PLM,gain,p,Om,W1,Wdiffs2,Wmid,Wmean,use_expectations_data,N);
    resnorm0 = sum(res0.^2);
    
    objh = @(alph) obj_GMM_LOMgain_univariate_mean(alph,x,fegrid_fine,param,gx,hx,eta,eN,vN,T,ndrop,PLM,gain,p,Om,W1,Wdiffs2,Wmid,Wmean,use_expectations_data,N);
    tic
    [alph_opt,resnorm,residual,flag] = lsqnonlin(objh,alph0,lb,ub,options);
    toc
    [res1, Om1, FE, Om_n, explode_t, negk_t] = obj_GMM_LOMgain_univariate_mean(alph_opt,x,fegrid_fine,param,gx,hx,eta,eN,vN,T,ndrop,PLM,gain,p,Om,W1,Wdiffs2,Wmid,Wmean,use_expectations_data,N);
    
    
    nancount = sum(sum(isnan(Om_n)));
    nanpercent = nancount/numel(Om_n);
    % treat the single outcomes as the mean outcomes
    alph_opt_mean = alph_opt;
    Om1mean = Om1;
    
    alpha_hat_boot(:,bb) = alph_opt; % NEW
end
