% command_sigmas.m
% check whether for different sigmas, the model can generate moments like
% those in the data
% 21 August 2020

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


skip = 1;
[fs, lw] = plot_configs;
datestr(now)

%% Compute weighting matrix and initialize alpha
% filename ='acf_data_11_Jun_2020'; % real data
% filename ='acf_data_21_Jul_2020'; % real data with SPF expectation in it
filename = 'acf_data_23_Aug_2020'; % real data with SPF expectation in it but qoq annualized inflation rates and expectations
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
% filename = 'acf_sim_univariate_data_10_Aug_2020'; % simulated data, nfe=5, fe=(-2,2), alph_true = 4*(0.05; 0.025; 0; 0.025; 0.05); referring to point (R,b) in my notes.
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

%% Generate moments for a single simulation of the true model


%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
alph_true = 4*[0.2,0.1,0,0.1,0.2]';
% alph_true = [1.0000    0.5000         0    0.5000    1.0000];
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
sig_r = 0.01;
sig_i = 0.01;
sig_u = 0.5;

eta = eye(nx).*[sig_r, sig_i, sig_u]'

param.psi_x = 0.3;
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

rng(2) %rng(0)
e = randn(nx,T+ndrop); % turned monpol shocks on in smat.m to avoid stochastic singularity!
v = 0*randn(ny+1,T+ndrop); % measurement error on the observables

[x0, y0, k0, phi0, FA0, FB0, FEt_10, diff0] = sim_learnLH_clean_approx_univariate(alph_true,x,param,gx,hx,eta, PLM, gain, T+ndrop,ndrop,e,v,knowTR,mpshock);

% This is the only new thing: add one-step ahead expectation to data (with iid shocks, it's just pibar, see Materials38, Point 3)
if use_expectations_data==1
    y_data = [y0(:,1:end-1); squeeze(phi0(1,1,1:end-1))'];
    % annualize inflation, interest rate and inflation expectations like in get_data.m
    y_data([1,3,4], :) = ((y_data([1,3,4], :)/100+1).^4 -1)*100;

elseif use_expectations_data==0
    % Don't use expectations
    y_data = y0(:,1:end-1);
    y_data([1,3], :) = ((y_data([1,3], :)/100+1).^4 -1)*100;
end
[nobs,T] = size(y_data)



% BK filter
K=12;
ystar_data = nan(nobs,T-2*K);
for i=1:nobs
    ystar_data(i,:) = BKfilter(y_data(i,:)');
end
lost_periods_BK = 2*K;


% choose your favorite filtered data
filt_data = ystar_data;
lost_periods = lost_periods_BK;


% compute moments, Om, as the autocovariances of the data for lags
% 0,1,...,K
K=4;%4
% Take the initial data, estimate a VAR
max_lags = 16;
[AIC,BIC,HQ] = aic_bic_hq(filt_data',max_lags);
p =min([AIC,BIC,HQ]);
% A is the impact matrix, identified via Cholesky, B is the beta_ols, res are
% the residuals, sigma is the estimated VC matrix.
[~,B,res,sigma] = sr_var(filt_data', p);

% return

% % ridge regression VAR
[B,~,sigma] = rf_var_ridge(filt_data', p, 0.001); 
% return

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
Om_model = vec(Gamj);


%% Plot covariogram of real data against that generated by the fake true model

% the appendix of each FIGNAME:
figspecs = ['N_', num2str(N),'_nfe_', num2str(nfe), ...
    '_gridspacing_', gridspacing, '_Wdiffs2_', num2str(Wdiffs2),'_Wmid_', num2str(Wmid), '_', cross_section, '_', this_code, '_', nowstr];

% Covariogram
Gamj_data = reshape(Om,nobs,nobs,K+1);
Gamj_model = reshape(Om_model,nobs,nobs,K+1);
% cvgram_data = zeros(nobs,K+1,nobs);
% cvgram0 = zeros(nobs,K+1,nobs);

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
        h = plot(0:K,squeeze(Gamj_data(i,j,:)), 'linewidth', lw);
        h0 = plot(0:K,squeeze(Gamj_model(i,j,:)), 'linewidth', lw);
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
    lh = legend([h,h0],{'Data', 'Model'},'interpreter', 'latex','Position',[0.45 -0.05 0.1 0.2], 'NumColumns',3, 'Box', 'off');
elseif contains(current_dir, 'gsfs0') % sirius server
    lh = legend([h,h0],{'Data', 'Model'},'interpreter', 'latex','Position',[0.45 -0.05 0.1 0.2], 'Box', 'off');
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