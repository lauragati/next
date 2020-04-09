% materials24
% 2 April 2020

% Goals:
% 1.) Figure out a way to guess a sequence of sthg and make it compatible
% with the model via optimization (target crit will be a special case)
% 2.) Implement 1.) using fsolve instead of fmincon
% 3.) Try the value function iteration to get at the optimal sequence
% 4.) Attempt Peter's idea of approximating the reaction function

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