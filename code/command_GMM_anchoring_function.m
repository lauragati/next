% command_GMM_anchoring_function.m
% a GMM implementation of estimating the anchoring function g
% created in materials22.m and replaces that file.
% 26 March 2020

clearvars
close all
clc

% Add all the relevant paths and grab the codename
this_code = mfilename;
[current_dir, basepath, BC_researchpath,toolpath,export_figpath,figpath,tablepath,datapath] = add_paths;
todays_date = strrep(datestr(today), '-','_');

% Variable stuff ---
print_figs        = 0;
stop_before_plots = 0;
skip_old_plots    = 0;
output_table = print_figs;

skip_old_stuff = 1;

%%

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
filt_data = ystar_data;

% 4.) Get autocovariances from 0 to lag K (override K)
K=4;
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
q=16; % blocksize 
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


% 6.) Do fmincon
%Optimization Parameters
options = optimset('fmincon');
options = optimset(options, 'TolFun', 1e-9, 'display', 'iter');

varp0 = [10,10]; % try estimating d and c
ub = [100,100];
lb = [0.1,0.1];
% %Compute the objective function one time with some values
% loss = objective_ACstructure(varp0,set,eN,T,N,burnin,PLM,gain, vec(acf_data),W1);

tic
%Declare a function handle for optimization problem
objh = @(varp) objective_ACstructure(varp,set,eN,T,N,burnin,PLM,gain, vec(acf_data));
[par_opt, loss_opt] = fmincon(objh, varp0, [],[],[],[],lb,ub,[],options);
% fmincon(FUN,X0,A,B,Aeq,Beq,LB,UB,NONLCON,OPTIONS)
toc

par_opt

loss_opt

% 7). Simulate the model using the optimal d:
param.d = par_opt(1);
param.c = par_opt(2);
[y,k] = fun_sim_anchoring(param,T,N, burnin,eN,PLM,gain);

if size(filt_data,2) == size(c_data,2)
    create_plot((1:numel(k)),1./k','k^{-1}',[this_code, '_gain_par_opt_HP', '_', todays_date],print_figs,'Inverse gain given estimated params')
elseif size(filt_data,2) == size(v_data,2)
    create_plot((1:numel(k)),1./k','k^{-1}',[this_code, '_gain_par_opt_Hamilton', '_', todays_date],print_figs,'Inverse gain given estimated params')
elseif size(filt_data,2) == size(ystar_data,2)
    create_plot((1:numel(k)),1./k','k^{-1}',[this_code, '_gain_par_opt_BK', '_', todays_date],print_figs,'Inverse gain given estimated params')
end