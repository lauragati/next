% A playground for value function iteration and policy iteration

clearvars
close all
clc

% Add all the relevant paths and grab the codename
this_code = mfilename;
[current_dir, basepath, BC_researchpath,toolpath,export_figpath,figpath,tablepath,datapath] = add_paths;
todays_date = strrep(datestr(today), '-','_');

skip = 1;


%% Peter's VFI example here, I'm also basing this on the Collard notes
if skip==0
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
    
end

%% Collard notes VFI
sig  = 1.5;
delt = 0.1;
bet  = 0.95;
alph = 0.3;
ngrid = 1000; % Collard's nbk
crit = 1;
epsi = 1e-6;

ks = ((1-bet*(1-delt))/(alph*bet))^(1/(alph-1));
dev = 0.9;
kmin = (1-dev)*ks;
kmax = (1+dev)*ks;
kgrid = linspace(kmin,kmax,ngrid)'; %% need to pay attention such that the dimensions of
% the grid, v, u, and c are the same (or at least compatible) for the
% maximzation, that's why I'm transposing here.
v = zeros(ngrid,1);
dr = zeros(ngrid,1); % decision rule = will store indices of optimal k' on the grid
vnew = zeros(ngrid,1);

iter = 0;
tic
while crit>epsi
iter = iter+1;
    for i=1:ngrid
        % select today's capital, k
        k = kgrid(i);
        % calculate the maximum index on the grid for k' such that
        % consumption is still positive
        imax = find(k^alph + (1-delt)*k - kgrid < 0,1) -1;
        if isempty(imax) % if it's empty, then take all the gridpoints
            imax = ngrid;
        end
        % now, given permissible values of k', calculate c...
        c = k^alph + (1-delt)*k - kgrid(1:imax);
        % ... and corresponding utility
        u = (c.^(1-sig)-1)/(1-sig);
        % update value function and choose index of optimal k' too
        [vnew(i), dr(i)] = max(u + bet*v(1:imax));        
    end
crit = max(abs(vnew-v));
v = vnew;
end
toc

% Final solution
kp = kgrid(dr);
c  = kgrid.^alph + (1-delt)*kgrid - kp;


figure
plot(kgrid,v)
title('Collard''s Fig. 7.1 Value function against k_t')

figure
subplot(1,2,1)
plot(kgrid,kp)
title('k_{t+1} as a function of k_t')
subplot(1,2,2)
plot(kgrid,c)
title('c_{t} as a function of k_t')
sgtitle('Collard''s Fig. 7.2 Decision rules')

%% Let's try to replicate Collard's "parametric dynamic programming example"
