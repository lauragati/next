% fun_GMM_LOMgain
% a function version of command_GMM_LOMgain so you can run it sequentially
% with different starting points
% 22 June 2020
function [alph_opt, resnorm, residual, flag] = fun_GMM_LOMgain(acf_outputs, nk1, nfe, k1min, k1max, femin, femax, alph0)



%% Compute weighting matrix and initialize alpha

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

% Specify info assumption on the Taylor rule and not to include a monpol
% shock - both set to 1 in the objective function
% knowTR =1  
% mpshock=1
%%%%%%%%%%%%%%%%%%%

% Size of cross-section
% we're not doing a whole cross-section here
ndrop = 5; % 0-50

% gen all the N sequences of shocks at once.
rng(0)
e = randn(ne,T+ndrop); % turned monpol shocks on in smat.m to avoid stochastic singularity!

% Create grids
% grids for k^(-1)_{t-1} and f_{t|t-1}
k1grid = linspace(k1min,k1max,nk1); 
fegrid = linspace(femin,femax,nfe); % for alph0, fe is between (-2.6278,3.5811). N=100 AR(1)-simulations of the model yield an average fe in (-0.2946,0.2809)

% map to ndim_simplex
x = cell(2,1);
x{1} = k1grid;
x{2} = fegrid;

% Let's plot the approximated evolution of the gain on a finer sample
ng_fine = 100;
k1grid_fine = linspace(k1min,k1max,ng_fine);
fegrid_fine = linspace(femin,femax,ng_fine);
[xxgrid_fine, yygrid_fine] = meshgrid(k1grid_fine,fegrid_fine);


%% GMM

% Fmincon
%Optimization Parameters
options = optimset('lsqnonlin');
options = optimset(options, 'TolFun', 1e-10, 'display', 'none'); 
options.MaxFunEvals = 15000;
options.MaxIter = 900;
options.UseParallel = 0; % 2/3 of the time


% let's keep these bounds
ub = k1max*ones(size(alph0));
lb = k1min*ones(size(alph0));

% %Compute the objective function one time with some values
% let's weight the prior...
Wprior=0; 

% tic
%Declare a function handle for optimization problem
objh = @(alph) obj_GMM_LOMgain(alph,x,xxgrid_fine,yygrid_fine,param,gx,hx,eta,e,T,ndrop,PLM,gain,p,Om,W1, alph0, Wprior);
[alph_opt,resnorm,residual,flag] = lsqnonlin(objh,alph0,lb,ub,options);
% toc




end

