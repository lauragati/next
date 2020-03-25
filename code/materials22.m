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

%2.) Set parameters for simulating anchoring model
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

% Size of cross-section
N = 1; % we're not doing a whole cross-section here
burnin = 100; % 100

% gen all the N sequences of shocks at once.
rng(0)
eN = randn(ne,T+burnin,N);


% 3.) Filter both using the same filter, at this point only the data
% 3.1) HP filter
g_data = nan(size(y_data));
c_data = nan(size(y_data));
for i=1:ny
    [g_data(i,:),c_data(i,:)] = HPfilter(y_data(i,:)');
end

% 3.2) Hamilton filter
h=8;
v_data = nan(ny,T-4-h+1);
for i=1:ny
    [v_data(i,:)] = Hamiltonfilter(y_data(i,:)');
end


% 3.3) BK filter
K=12;
ystar_data = nan(ny,T-2*K);
for i=1:ny
    ystar_data(i,:) = BKfilter(y_data(i,:)');
end

% % Plot filtered inflation
% create_plot(1:numel(c_data(1,:)),c_data(1,:),'HP-filtered cycle',[this_code, '_HP'],print_figs,'PCE Inflation')
% create_plot(1:numel(v_data(1,:)),v_data(1,:),'Hamilton-filtered cycle',[this_code, '_Hamilton'],print_figs,'PCE Inflation')
% create_plot(1:numel(ystar_data(1,:)),ystar_data(1,:),'BK-filtered cycle',[this_code, '_BK'],print_figs,'PCE Inflation')


% choose your favorite filtered data
filt_data = v_data;

% 4.) Get autocovariances from 0 to lag K (override K)
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

acf_data = zeros(ny,K+1);
for i=1:ny
    for j=1:K+1
        k=j-1;
        VC_data = cov(filt_data(i,k+1:end), filt_data(i,1:end-k));
        acf_data(i,j) =VC_data(2,1);
    end
end

% 5.) Get bootstrapped variance of data moments
nboot =10000;
q=8; % blocksize 
resamples = bootstrap_resample(filt_data, nboot, q);

acf_data_boot = zeros(ny,K+1,nboot);
for n = 1:nboot
    for i=1:ny
        for j=1:K+1
            k=j-1;
            VC_data_boot = cov(resamples(i,k+1:end,nboot), resamples(i,1:end-k,nboot));
            acf_data_boot(i,j,n) =VC_data_boot(2,1);
        end
    end
end
% I'm surprised: all of this in 5) takes less than 4 sec.
W = zeros(1,ny*(K+1));
index=0;
for i=1:ny
     for j=1:K+1
         index=index+1;
         W(index) = cov(squeeze(acf_data_boot(i,j,:)));
     end
end
W1 = W.^(-1);
W1 = diag(W1);

return
% 6.) Do fmincon
%Optimization Parameters
options = optimset('fmincon');
options = optimset(options, 'TolFun', 1e-9, 'display', 'iter');

varp0 = 10; % right now I'm only estimating d
ub = 100;
lb = 0.1;
% %Compute the objective function one time with some values
% loss = objective_ACstructure(varp0,set,eN,T,N,burnin,PLM,gain, vec(acf_data),W1);

tic
%Declare a function handle for optimization problem
objh = @(varp) objective_ACstructure(varp,set,eN,T,N,burnin,PLM,gain, vec(acf_data));
[par_opt, loss_opt] = fmincon(objh, varp0, [],[],[],[],lb,ub,[],options);
% fmincon(FUN,X0,A,B,Aeq,Beq,LB,UB,NONLCON,OPTIONS)
toc

par_opt

% 7). Simulate the model using the optimal d:
param.d = par_opt;
[y,k] = fun_sim_anchoring(param,T,N, burnin,eN,PLM,gain);

if size(filt_data,2) == size(c_data,2)
create_plot((1:numel(k)),1./k','k^{-1}',[this_code, '_gain_dhat_HP'],print_figs,'Inverse gain given estimated d')
elseif size(filt_data,2) == size(v_data,2)
    create_plot((1:numel(k)),1./k','k^{-1}',[this_code, '_gain_dhat_Hamilton'],print_figs,'Inverse gain given estimated d')
elseif size(filt_data,2) == size(ystar_data,2)
    create_plot((1:numel(k)),1./k','k^{-1}',[this_code, '_gain_dhat_BK'],print_figs,'Inverse gain given estimated d')
end

%% 2.) Numerical optimal plan
clearvars
clc
[param, set, param_names, param_values_str, param_titles] = parameters_next;

% Structure will be
% 1.) Generate a guess exog interest rate sequence
T = 100; H = 50; burnin = 0; N=1; ne=3;
i_seq0 = gen_AR1(T+H,0.9,1);

% and gen shocks
rng(0)
eN = randn(ne,T+H+burnin,N);

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

% 2.) Solve for optimal interest rate path using a target-criterion-based
% loss.
% I'll call this most optimal of plans the Ramsey plan. 

%Optimization Parameters
options = optimset('fmincon');
options = optimset(options, 'TolFun', 1e-9, 'display', 'iter');

ub = 40*ones(T,1);
lb = 0.001*ones(T,1);
% %Compute the objective function one time with some values
loss = objective_target_criterion(i_seq0,param,eN,T,N,burnin,PLM,gain);

tic
%Declare a function handle for optimization problem
objh = @(varp) objective_target_criterion(varp,param,eN,T,N,burnin,PLM,gain);
[i_ramsey, loss_opt] = fmincon(objh, i_seq0, [],[],[],[],lb,ub,[],options);
% fmincon(FUN,X0,A,B,Aeq,Beq,LB,UB,NONLCON,OPTIONS)
toc

i_ramsey

% 3.) Simulate model given the Ramsey plan for i -> obtain Ramsey plans for
% x and pi
[pi_ramsey,x_ramsey,k] = fun_sim_anchoring_given_i(param,T+H,N,burnin,eN,PLM,gain,i_ramsey);


% see (plot) deviations between the Ramsey plan and the plans obtained when
% using a Taylor rule
[y_TR,k_TR] = fun_sim_anchoring(param,T+H,N, burnin,eN,PLM,gain);

% 4.) Can even do search over TR parameters to see if they can implement
% the Ramsey plan.