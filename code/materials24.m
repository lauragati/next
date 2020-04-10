% materials24
% 2 April 2020

% Goals:
% 1.) Figure out a way to guess a sequence of sthg and make it compatible
% with the model via optimization (target crit will be a special case)
% 2.) Implement 1.) using fsolve instead of fmincon
% 3.) Try the value function iteration to get at the optimal sequence
% 4.) Attempt Peter's idea of approximating the reaction function
% 5.) create sim_learnLH_clean.m, a cleaned-up but identical version of
% sim_learnLH.m

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

skip = 1;
if skip==0
    %% 1) Simulate given a sequence - optimize over that sequence to satisfy model
    command_sim_given_seq
    % this needs to be corrected, there is an fsolve way to do it conceptually
    % better
    % 2)
    command_sim_given_seq_fsolve
    
    %%  Value function iteration to solve for optimal i-sequence  ain't workin yet
    command_valfun_iter
end


%% 3.) Work thru Peter's VFI example here, I'm also basing this on the Collard notes
play_value_function_iter

clc
% params
tol=0.01; maxiter =200; dif = tol+1000; iter=0;
bet=0.99; alph = 0.3;
ks = (1/(alph*bet))^(1/(alph-1));

dev=0.9;
kmin = (1-dev)*ks;
kmax = (1+dev)*ks;
ngrid = 1000; % number of data points in the grid (nbk in Collard's notation)
dk = (kmax-kmin)/(ngrid-1); % implied increment
kgrid = linspace(kmin,kmax,ngrid);

v = zeros(1, ngrid);
vnew = zeros(1,ngrid);
jstar = zeros(1,ngrid);
tic
while dif > tol && iter < maxiter
    for i=1:ngrid
        % The idea is to do the max for each k_t(i)
        k = kgrid(i);
        % First, for this value of k_t, k_{t+1} must be less than k_t^alph
        ub = k^alph;
        % find the lowest index at which k_{t+1} exceeds k_t^alph
        imax = find(kgrid > ub, 1 );
        if isempty(imax) % if it's empty, then take all the gridpoints
            imax = ngrid;
        end
        % 2nd, cons and utility as a function of this k_t, for all for all values of k_{t+1}
        c = k^alph - kgrid(1:imax); % c(i,:) for all values of k_{t+1}
        u = log(c); % u(i,:) for all values of k_{t+1}
        [vnew(i), jstar(i)] = max(u + bet*v(1:imax));
    end
    dif = max(abs((vnew-v)));
    v = vnew;
    iter = iter+1;
end
toc

% Final sol
% 1.) k_{t+1} as a function of k_t
kp = kgrid(jstar);
% 2.) c as function of k_t, k_{t+1}
c = kgrid.^alph - kp;


figure
plot(kgrid,c)
title('Consumption as a function of k_t')

figure
plot(kgrid,kp)
title('k_{t+1} as a function of k_t')

%% 4.) approximating the reaction function
if skip==0
    command_approx_reaction
end

%% 5.) Clean up sim_learnLH.m - this here copies % command_IRFs_anchoring
% Do IRFs for anchoring true-baseline learning models
% adapted from command_IRFs_many_learning.m
% 25 Jan 2020

clearvars
close all
clc

date_today = strrep(datestr(today, 'yyyy-mm-dd'), '-','_');

% % Add all the relevant paths and grab the codename
this_code = mfilename;
[current_dir, basepath, BC_researchpath,toolpath,export_figpath,figpath,tablepath,datapath] = add_paths;

% Variable stuff ---
print_figs        = 0;
stop_before_plots = 0;
skip_old_plots    = 0;
output_table = print_figs;

plot_IRFs=1;
plot_simulated_sequence = 1;
plot_gains=0;
plot_gain_IRF = 0;
plot_IRFs_anch = 0;
skip_old_stuff = 1;

% Parameters
tic
[param, set, param_names, param_values_str, param_titles] = parameters_next;


psi_x = param.psi_x;
psi_pi = param.psi_pi;
thetbar =param.thetbar;
thettilde = param.thettilde;

sig_r = param.sig_r;
sig_i = param.sig_i;
sig_u = param.sig_u;
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


% Model selection and informational assumption

% Model selection
%%%%%%%%%%%%%%%%%%%
PLM = constant_only;
gain = dgain;
%%%%%%%%%%%%%%%%%%%
% the smooth criterion only works in this combo:
%%%%%%%%%%%%%%%%%%%
PLM = constant_only_pi_only;
gain = again_critsmooth;
%%%%%%%%%%%%%%%%%%%

T = 400 % 400
% Size of cross-section
N = 100;%100 500
burnin = 100;
dt_vals = 25; %25 time of imposing innovation 345
h = 10; % h-period IRFs

% RE model
[fyn, fxn, fypn, fxpn] = model_NK(param);
[gx,hx]=gx_hx_alt(fyn,fxn,fypn,fxpn);
[ny, nx] = size(gx);
[Aa, Ab, As] = matrices_A_13_true_baseline(param, hx);
SIG = eye(nx).*[sig_r, sig_i, sig_u]';
eta = SIG; %just so you know

[PLM_name, gain_name, gain_title] = give_names(PLM, gain);


% Simulate models

% gen all the N sequences of shocks at once.
rng(0)
eN = randn(ne,T+burnin,N);

% Preallocate
nd = size(dt_vals,2);
d = 1; % the innovation, delta
GIR_Y_LH = zeros(ny,h,N,nd);
GIR_k = zeros(h,N,nd);
k = zeros(T,N);
ks= zeros(T,N);
k_clean = zeros(T,N);
k_smooth = zeros(T,N);
anch = zeros(1,N); % indicator for whether that sequence was anchored when the shock hit
g_pi = zeros(T,1);
g_pibar = zeros(T,1);
% warning off
for s=2  %2->zoom in on monetary policy shock
    x0 = zeros(1,nx);
    x0(s) = d;
    
    for n=1:N
        % Sequence of innovations
        e = squeeze(eN(:,:,n));
        
        % Unshocked
        %         dbstop if warning
        % RE
        [x_RE, y_RE] = sim_model(gx,hx,SIG,T,burnin,e);
        % Learning
        [x_LH, y_LH, ~, ~, ~, ~, ~, ~, diff(:,n),~, k(:,n)] = sim_learnLH(gx,hx,SIG,T+burnin,burnin,e, Aa, Ab, As, param, PLM, gain);
        [x_clean, y_clean, k_clean(:,n), ~, ~, ~, diff_clean(:,n)] = sim_learnLH_clean(param,gx,hx,SIG, Aa, Ab, As,PLM, gain, T+burnin,burnin,e);
        [x_smooth, y_smooth, k_smooth(:,n), pibar(:,n), FA, FB, g_pi, g_pibar, diff_smooth(:,n)] = sim_learnLH_clean_smooth(param,gx,hx,SIG, Aa, Ab, As,T+burnin,burnin,e);
        % Shocked
        % RE
        % make RE shock the same scale as learning:
        x0RE = (SIG*x0')';
        [IR, iry, irx]=ir(gx,hx,x0RE,h);
        iry = iry';
        RE_fcsts = gx*hx*irx';
        
        % Learning
        for t=1:nd
            dt = dt_vals(t);
            
            % Shocked
            [~, ys_LH, ~, ~, ~, ~, ~, ~, ~,~, ~] = sim_learnLH(gx,hx,SIG,T+burnin,burnin,e, Aa, Ab, As, param, PLM, gain, dt, x0);
            [~, ys_clean, ~, ~, ~, ~, ~] = sim_learnLH_clean(param,gx,hx,SIG, Aa, Ab, As,PLM, gain, T+burnin,burnin,e, dt, x0);
            [~, ys_smooth, ~, ~, ~, ~, ~] = sim_learnLH_clean_smooth(param,gx,hx,SIG, Aa, Ab, As,T+burnin,burnin,e,dt,x0);
        end
        
    end
end

disp(['Max difference in states = ', num2str(max(max(abs(x_LH-x_clean))))])
disp(['Max difference in jumps = ', num2str(max(max(abs(y_LH-y_clean))))])
disp(['Max difference in k = ', num2str(max(max(abs(k-k_clean))))])
disp(['Max difference in convergence diffs = ', num2str(max(max(abs(diff-diff_clean))))])
disp(['Max difference in shocked jumps = ', num2str(max(max(abs(ys_LH-ys_clean))))])
disp('Diffs for sim_learnLH_clean_smooth.m')
disp(['Max difference in states = ', num2str(max(max(abs(x_LH-x_smooth))))])
disp(['Max difference in jumps = ', num2str(max(max(abs(y_LH-y_smooth))))])
disp(['Max difference in k = ', num2str(max(max(abs(k-k_smooth))))])
disp(['->Max difference in k_clean, k_smooth = ', num2str(max(max(abs(k_clean-k_smooth))))])
disp(['Max difference in convergence diffs = ', num2str(max(max(abs(diff-diff_smooth))))])
disp(['Max difference in shocked jumps = ', num2str(max(max(abs(ys_LH-ys_smooth))))])


figure
plot(mean(diff,2))
title('Convergence - mean(diff) across N')

figure
plot(diff(:,end))
title('Convergence  - diff for N=end')

figure
plot(mean(pibar,2))
title('mean(\bar{\pi})')

figure
plot(pibar(:,end))
title('Pibar  - diff for N=end')


