% fun_GMM_LOMgain_univariate
% command_GMM_LOMgain_univariate
% Same as command_GMM_LOMgain_univariate, except as a function so you can
% start at different starting values
% 22 June 2020

function [alph_opt, resnorm, residual, Om_opt, flag] = fun_GMM_LOMgain_univariate(acf_outputs, nfe, k1min,k1max,femin, femax, alph0)

%% Compute weighting matrix and initialize alpha

%%%%%%%%%%%%%%%%%%%
% Grid

% upper and lower bounds for estimated coefficients
ub = k1max*ones(nfe,1); %1
lb = k1min*ones(nfe,1); %0
% weights on additional moments
Wprior=0;%0
Wconvexity=1000;%1000
Wmean=100;%0


%Optimization Parameters
options = optimoptions('lsqnonlin');
options = optimoptions(options, 'display', 'none');
options.TolFun= 1e-9;
% options.OptimalityTolerance = 1e-9; % this is the guy you can access in optimoptions, not in optimset. It pertains to first order optimality measure.
options.MaxFunEvals = 1000;
% options.MaxIter = 1200;
options.TolX = 1e-9;
options.UseParallel = 1; % 2/3 of the time
%%%%%%%%%%%%%%%%%%%



Om = acf_outputs{1}; % this is the moments vector
Om_boot = acf_outputs{2}; % moments vectors in bootstrapped samples
ny = acf_outputs{3};
p = acf_outputs{4};
K = acf_outputs{5};
filt_data = acf_outputs{6};
lost_periods = acf_outputs{7};
T = length(filt_data);
T = T+lost_periods; % to make up for the loss due to filtering


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
param = parameters_next;
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


% Size of cross-section
% we're not doing a whole cross-section here
ndrop = 5; % 0-50

% gen all the N sequences of shocks at once.
rng(0)
e = randn(ne,T+ndrop); % turned monpol shocks on in smat.m to avoid stochastic singularity!

fegrid = linspace(femin,femax,nfe); % for alph0, fe is between (-2.6278,3.5811).
% map to ndim_simplex
x = cell(1,1);
x{1} = fegrid;

% Let's plot the initial approximated evolution of the gain on a finer sample
ng_fine = 100;
fegrid_fine = linspace(femin,femax,ng_fine);


%% GMM

%Declare a function handle for optimization problem
objh = @(alph) obj_GMM_LOMgain_univariate(alph,x,fegrid_fine,param,gx,hx,eta,e,T,ndrop,PLM,gain,p,Om,W1,Wconvexity,Wmean);
[alph_opt,resnorm,residual,flag] = lsqnonlin(objh,alph0,lb,ub,options);


% Plot ACFs at start and end (Om0 and Om1 are the model-implied moments, initial and optimal)
[~, Om_opt] = obj_GMM_LOMgain_univariate(alph_opt,x,fegrid_fine,param,gx,hx,eta,e,T,ndrop,PLM,gain,p,Om,W1,Wconvexity,Wmean);

end


