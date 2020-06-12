% materials33_server.m
%% Compute weighting matrix and do GMM
filename ='acf_data_11_Jun_2020';
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

% weighting matrix for GMM
% note: 2nd argument of var(X,W,DIM) signifies whether to normalize by N-1
% or N. 0 -> N-1, 1 -> N. The third argument, DIM, says along which
% dimension to take the variance.
W = diag(var(Om_boot,0,2));
W1 = W^(-1);

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
PLM = constant_only;
gain = again_critsmooth;
%%%%%%%%%%%%%%%%%%%

% Size of cross-section
% we're not doing a whole cross-section here
%%%%%%%%%%%
ndrop = 50; % 100
%%%%%%%%%%%

% gen all the N sequences of shocks at once.
rng(0)
e = randn(ne,T+ndrop); % turned monpol shocks on in smat.m to avoid stochastic singularity!

% Fmincon
%Optimization Parameters
options = optimset('fmincon');
options = optimset(options, 'TolFun', 1e-9, 'display', 'iter'); 
% options.MaxFunEvals = 10000;
options.UseParallel = 1; % 2/3 of the time

% Do an initial approx of the anchoring function to initialize the coeffs
ng = 10;
% grids for k^(-1)_{t-1} and f_{t|t-1}
k1grid = linspace(0.001,param.gbar,ng);
fegrid = linspace(-5,5,ng);
% values for k^{-1}_t for the grid
k = zeros(ng,ng);
for i=1:ng
    for j=1:ng
        k(i,j) = fk_smooth_pi_only(param,fegrid(j), 1./k1grid(i));
    end
end
% map to ndim_simplex
x = cell(2,1);
x{1} = k1grid;
x{2} = fegrid;
[xxgrid, yygrid] = meshgrid(k1grid,fegrid);
% % check quikcly whether it looks ok
kmesh = fk_smooth_pi_only(param,yygrid,1./xxgrid);
disp(['max(kmesh-k'') = ', num2str(max(max(abs(kmesh - k'))))])
k1 = 1./kmesh;

alph0 = ndim_simplex(x,[xxgrid(:)';yygrid(:)'],k1);

% Let's plot the approximated evolution of the gain on a finer sample
ng_fine = 100;
k1grid_fine = linspace(0.001,param.gbar,ng_fine);
fegrid_fine = linspace(-5,5,ng_fine);
[xxgrid_fine, yygrid_fine] = meshgrid(k1grid_fine,fegrid_fine);

k10 = ndim_simplex_eval(x,[xxgrid_fine(:)';yygrid_fine(:)'],alph0);



ub = alph0+0.2;
lb = zeros(size(alph0));
% %Compute the objective function one time with some values
loss = obj_GMM_LOMgain(alph0,x,param,gx,hx,eta,e,T,ndrop,PLM,gain,Om,W1);

tic
%Declare a function handle for optimization problem
objh = @(alph) obj_GMM_LOMgain(alph,x,param,gx,hx,eta,e,T,ndrop,PLM,gain,Om,W1);
[alph_opt, loss_opt] = fmincon(objh, alph0+0*rand(size(alph0)), [],[],[],[],lb,ub,[],options);
% fmincon(FUN,X0,A,B,Aeq,Beq,LB,UB,NONLCON,OPTIONS)
toc

% Let's add the final output to the finer sample
k1_opt = ndim_simplex_eval(x,[xxgrid_fine(:)';yygrid_fine(:)'],alph_opt);

filename = ['estim_outputs_server_nburn_', num2str(ndrop), '_date_', todays_date];
estim_outputs = {xxgrid_fine,yygrid_fine,ng_fine,k10,k1_opt,filename};
save([filename,'.mat'], 'estim_outputs')
