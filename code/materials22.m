% materials22
% 19 March 2020

% Goals:
% 1.) Do initial GMM estimation of the anchoring function in midsimple
% model - think I can take some shortcuts from general codes
% 2.) Do the numerical implementation of the target criterion in midsimple
% model - same shortcuts apply?

clearvars
close all
clc

% Add all the relevant paths and grab the codename
this_code = mfilename;
[current_dir, basepath, BC_researchpath,toolpath,export_figpath,figpath,tablepath,datapath] = add_paths;

% Variable stuff ---
print_figs        = 0;
stop_before_plots = 0;
skip_old_plots    = 0;
output_table = print_figs;

skip_old_stuff = 1;

%% 1.) GMM
% 1.) Get real data
% Suppose we have real data (in logs)
y_data = get_data;
[ny,T] = size(y_data)

%2.) Simulate anchoring model
[param, set, param_names, param_values_str, param_titles] = parameters_next;
ne = 3;

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

% T = 100 % 400
% Size of cross-section
N = 1; % we're not doing a whole cross-section here
burnin = 100; % 100

% gen all the N sequences of shocks at once.
rng(0)
eN = randn(ne,T+burnin,N);

% So the simulated data are:
y = fun_sim_anchoring(param,T,N, burnin,eN,PLM,gain);
ny = size(y,1);


% 2.) Filter both using the same filter
% 2.1) HP filter
g = nan(size(y));
c = nan(size(y));
g_data = nan(size(y));
c_data = nan(size(y));
for i=1:ny
    [g(i,:),c(i,:)] = HPfilter(y(i,:)');
    [g_data(i,:),c_data(i,:)] = HPfilter(y_data(i,:)');
end

% 2.2) Hamilton filter
v = nan(size(y));
v_data = nan(size(y));
for i=1:ny
    [v(i,:)] = HPfilter(y(i,:)');
    [v_data(i,:)] = HPfilter(y_data(i,:)');
end

% 2.3) BK filter
K=12;
ystar = nan(ny,T-2*K-1);
ystar_data = nan(ny,T-2*K-1);
for i=1:ny
    ystar(i,:) = BKfilter(y(i,:)');
    ystar_data(i,:) = BKfilter(y_data(i,:)');
end

% choose your favorite filtered data
filt = ystar;
filt_data = ystar_data;

% 3.) Get autocovariances from 0 to lag K (override K)
K=4;
% A nice check because we know the autocov of an AR1
% ar1 = gen_AR1(T, 0.9, 1);
% cov(ar1(2:end), ar1(1:end-1))
% corr(ar1(2:end), ar1(1:end-1)) % note: corr = cov/(sqrt(var(x)*var(y)))
% corr(ar1(3:end), ar1(1:end-2))
% for i=1:K+1
%     k = i-1;
%     VC = cov(ar1(k+1:end), ar1(1:end-k));
%     acf(i) = VC(2,1);
% end

acf = zeros(ny,K+1);
acf_data = zeros(ny,K+1);
for i=1:ny
    for j=1:K+1
        k=j-1;
        VC = cov(filt(i,k+1:end), filt(i,1:end-k));
        acf(i,j) =VC(2,1);
        VC_data = cov(filt_data(i,k+1:end), filt_data(i,1:end-k));
        acf_data(i,j) =VC_data(2,1);
    end
end

% 4.) Do fmincon
%Optimization Parameters
options = optimset('fmincon');
options = optimset(options, 'TolFun', 1e-9, 'display', 'iter'); 

varp0 = 10; % right now I'm only estimating d
ub = 100; 
lb = 0.1;
%Compute the objective function one time with some values
loss = objective_ACstructure(varp0,set,eN,T,N,burnin,PLM,gain, vec(acf_data));

tic
%Declare a function handle for optimization problem
objh = @(varp) objective_ACstructure(varp,set,eN,T,N,burnin,PLM,gain, vec(acf_data));
[par_opt, loss_opt] = fmincon(objh, varp0, [],[],[],[],lb,ub,[],options);
% fmincon(FUN,X0,A,B,Aeq,Beq,LB,UB,NONLCON,OPTIONS)
toc

par_opt





