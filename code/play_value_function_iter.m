% A playground for value function iteration

clearvars
close all
clc

% Add all the relevant paths and grab the codename
this_code = mfilename;
[current_dir, basepath, BC_researchpath,toolpath,export_figpath,figpath,tablepath,datapath] = add_paths;
todays_date = strrep(datestr(today), '-','_');

print_figs=0;
skip = 1;



if skip==0
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
    figname = 'vfi_discrete_valuefun';
    
    if print_figs ==1
        disp(figname)
        cd(figpath)
        export_fig(figname)
        cd(current_dir)
        close
    end
    
    
    figure
    subplot(1,2,1)
    plot(kgrid,kp)
    title('k_{t+1} as a function of k_t')
    subplot(1,2,2)
    plot(kgrid,c)
    title('c_{t} as a function of k_t')
    sgtitle('Collard''s Fig. 7.2 Decision rules')
    % it works!
    
    figname = 'vfi_discrete_decisions';
    
    if print_figs ==1
        disp(figname)
        cd(figpath)
        export_fig(figname)
        cd(current_dir)
        close
    end
    
    %% Let's try to replicate Collard's "parametric dynamic programming example"
    sig  = 1.5;
    delt = 0.1;
    bet  = 0.95;
    alph = 0.3;
    
    ngrid=20; % Collard's nbk, he says he set it to 20, but maybe not...
    % p = 10; % order of Chebyshev polyonomials (0,1,...,10)
    crit=1;
    iter=1;
    epsi=1e-6;
    ks = ((1-bet*(1-delt))/(alph*bet))^(1/(alph-1));
    dev = 0.9;
    kmin = (1-dev)*ks;
    kmax = (1+dev)*ks;
    zk =  -cos((2*[1:ngrid]'-1)*pi/(2*ngrid));
    kgrid = (kmax-kmin)*(zk+1)/2 +kmin;
    % Initial guess for the approximation
    % Collard calls the coefficients theta, thus th0
    % Collard's logic is: initialize value function at utility function...
    v       = (((kgrid.^alph).^(1-sig)-1)/((1-sig)*(1-bet)));
    % ...calculate Chebyshev polynomials at the [-1,1] nodes...
    X       = chebychev(zk,ngrid);
    % ...and back out initial coefficients from the approximation identity V =thet*cheby
    th0     = X\v; % initial guess for params
    Tv      = zeros(ngrid,1);
    kp      = zeros(ngrid,1);
    
    
    %Optimization Parameters
    options = optimset('fmincon');
    options = optimset(options, 'TolFun', 1e-9, 'display', 'none');
    
    tic
    warning off
    datestr(now)
    while crit>epsi
        k0 = kgrid(1);
        for i=1:ngrid
            param = [alph, bet, delt, sig, kmin, kmax, ngrid, kgrid(i)];
            % kp(i) = some optimization routine that gives you optimal k'
            % Given current coefficients th0 and today's k, kgrid(i), starting
            % with a guess of k0, solve for k':
            objh = @(kp) tv(kp,param,th0);
            %         [kp(i)] = fmincon(objh, k0, [],[],[],[],kmin,kmax,[],options);
            [kp(i)] = fminunc(objh, k0,options);
            % fmincon(FUN,X0,A,B,Aeq,Beq,LB,UB,NONLCON,OPTIONS)
            % I wonder if we could get around doing fmincon...
            
            k0 = kp(i);
            %get back today's value function given the optimal kp(i) and the
            %T-map given current "estimated" coefficients theta0:
            Tv(i) = -tv(kp(i), param, th0);
        end
        % Fitting step: Given updated value function, back out coefficients
        thet =X\Tv;
        crit = max(abs(Tv-v));
        v = Tv;
        th0 = thet;
        iter=iter+1;
    end
    toc % indeed it takes 242 iterations (But it takes 60 sec. With fminunc only 12 sec!)
    
    % Compute consumption function
    c = kgrid.^alph+(1-delt)*kgrid-kp;
    
    figure
    plot(kgrid,v)
    title('Collard''s Fig. 7.7 Value function against k_t')
    
    figure
    subplot(1,2,1)
    plot(kgrid,kp)
    title('k_{t+1} as a function of k_t')
    subplot(1,2,2)
    plot(kgrid,c)
    title('c_{t} as a function of k_t')
    sgtitle('Collard''s Fig. 7.8 Decision rules')
    % it works!
    
    
    %% Let's try the parametric VFI on the optimal growth with Collard's parameterization, but on my own
    % let me use notation in Cai and Judd 2014
    clc
    
    sig  = 1.5;
    delt = 0.1;
    bet  = 0.95;
    alph = 0.3;
    param = [sig,delt,bet,alph];
    
    %Optimization Parameters
    % options = optimset('fmincon');
    options = optimset('fminunc');
    options = optimset(options, 'TolFun', 1e-16, 'display', 'none');
    
    m=20; % number of nodes (index i)
    n = 19; % 10 order of Chebyshev polyonomials (0,1,...,10) (index j)
    crit=1;
    iter=1;
    maxiter=300;%300
    epsi=1e-6; %-6
    ks = ((1-bet*(1-delt))/(alph*bet))^(1/(alph-1));
    dev = 0.9;
    kmin = (1-dev)*ks;
    kmax = (1+dev)*ks;
    zi =  -cos((2*(1:m)'-1)*pi/(2*m));
    xgrid = (zi+1).*(kmax-kmin)/2 +kmin;
    % cheby polynomials on this grid (are constant y'know)
    Tn = chebyshev(zi,n);
    % Initial nodes as a completely uneducated guess
    b0 = ones(n+1,1);
    v0 = vhat_cheby(b0,xgrid,n,kmin,kmax);
    % % initialize the value function instead as umax given kp=0
    % v0 = (((xgrid.^alph).^(1-sig)-1)/((1-sig)*(1-bet))); % -> explodes!
    % b0 =Tn\v0;
    v1 = zeros(size(v0));
    kp = zeros(m,1);
    
    b=b0;
    v=v0;
    
    % always start searching at the beginning of the grid
    k0 = xgrid(1);
    datestr(now)
    tic
    while crit > epsi && iter< maxiter
        % Step 1: maximization
        for j=1:m
            kj=xgrid(j);
            negv = @(kp) TV(kp,b,xgrid(j),n,kmin,kmax, param);
            
            %         [kp(j)] = fmincon(negv, k0, [],[],[],[],kmin,kmax,[],options);
            [kp(j,iter)] = fminunc(negv, k0, options); % noticeably faster
            % compute the value function at the maximizing kp
            v1(j) = -TV(kp(j,iter),b,kj,n,kmin,kmax, param);
            
            %         % let's also do it with Collard's methods
            %         negv_collard = @(kp) tv(kp,[alph, bet, delt, sig, kmin, kmax, n, kj],b);
            %         [kp_collard(j,iter)] = fminunc(negv_collard, k0, options); % these equal mine
            %         v1_collard(j) = -tv(kp(j,iter),[alph, bet, delt, sig, kmin, kmax, n, kj],b); % these equal mine too
            
            
            % maybe this makes fmincon search elsewhere? It does, it seems like
            % an educated guess. It does make a diff too.
            %         k0 = kp(j);
            %         disp(kp(j))
        end
        
        % Step 2: fitting
        %     b1 = compute_cheby_coeffs(v,xgrid,n,kmin,kmax); % don't actually need
        %     this!
        b1 = Tn\v1; % b1-b1_alt not 0.
        
        % Compute stopping criterion and update
        crit = max(abs(v1-v));
        if mod(iter,50)==0
            disp(['Concluding iteration = ', num2str(iter)])
            disp(['Stopping criterion = ', num2str(crit)])
        end
        iter=iter+1;
        %     b=b1;
        b=b1;
        v=v1;
    end
    toc
    
    % Compute consumption function
    c = xgrid.^alph+(1-delt).*xgrid-kp(:,end);
    
    
    figure
    plot(xgrid,v)
    title('Collard''s Fig. 7.7 Value function against k_t - Chebyshev polynomials')
    figname = 'vfi_cheby_valuefun';
    if print_figs ==1
        disp(figname)
        cd(figpath)
        export_fig(figname)
        cd(current_dir)
        close
    end
    
    figure
    subplot(1,2,1)
    plot(xgrid,kp(:,end)); hold on
    plot(xgrid,xgrid,'k--')
    title('k_{t+1} as a function of k_t')
    subplot(1,2,2)
    plot(xgrid,c)
    title('c_{t} as a function of k_t')
    sgtitle('Collard''s Fig. 7.8 Decision rules - Chebyshev polynomials')
    figname = 'vfi_cheby_decisions';
    
    if print_figs ==1
        disp(figname)
        cd(figpath)
        export_fig(figname)
        cd(current_dir)
        close
    end
    
    % it works, yeaaaaah! and takes <9 sec.
    
    %% some initial testing of functions
    zi =  -cos((2*[1:m]'-1)*pi/(2*m));
    kgrid = (kmax-kmin)*(zi+1)/2 +kmin;
    
    % 1.) Chebyshev polyonomials should equal
    T_me = chebyshev(zi,n);
    T_collard = chebychev(zi,n);
    T_me-T_collard % yep
    
    % 2.) vhat given coefficients should equal
    v0 = (((kgrid.^alph).^(1-sig)-1)/((1-sig)*(1-bet)));% take Collard's initial value function and coefficients
    b0 =T_me\v0;
    vhat_collard = value(ks,[kmin, kmax, n],b0);
    vhat_me = vhat_cheby(b0,ks,n,kmin,kmax);
    vhat_collard - vhat_me % yep
    
    % 3.) TV
    negv_collard = tv(ks,[alph, bet, delt, sig, kmin, kmax, n, ks],b0);
    negv_me = TV(ks,b0,ks,n,kmin,kmax, [sig,delt,bet,alph]);
    negv_collard-negv_me
    
    
    %
    % vhat1 = vhat_cheby(ones(1,n+1),xgrid,n,kmin,kmax);
    % % plot(vhat1) % it's very oscillatory, no idea if that's a good thing
    %
    % tvhat1 = TV(xgrid(1),b0,xgrid(1),n,kmin,kmax, param); % ok seems to work
    %
    % b1 = compute_cheby_coeffs(vhat1,xgrid,n,kmin,kmax); % ok seems to work
    
    %% let's try parametric VFI with a spline!
    
    clc
    
    sig  = 1.5;
    delt = 0.1;
    bet  = 0.95;
    alph = 0.3;
    param = [sig,delt,bet,alph];
    
    %Optimization Parameters
    % options = optimset('fmincon');
    options1 = optimset('fminunc');
    options1 = optimset(options1, 'TolFun', 1e-16, 'display', 'none');
    options2 = optimset('fsolve');
    options2 = optimset(options2, 'TolFun', 1e-16, 'display', 'none');
    
    m=20; % number of nodes (index i), m-1 intervals
    crit=1;
    iter=1;
    maxiter=600;%300
    epsi=1e-6; %-6
    ks = ((1-bet*(1-delt))/(alph*bet))^(1/(alph-1));
    dev = 0.9;
    kmin = (1-dev)*ks;
    kmax = (1+dev)*ks;
    % for the spline, uniformly placed nodes are fine
    xgrid = linspace(kmin,kmax,m)';
    % Initial nodes as a completely uneducated guess
    % coeffs0 = zeros(m-1,4); % [a,b,c,d] % this actually converges the fastest, but you always get the same thing
    % coeffs0 = ones(m-1,4);
    % coeffs0 = rand(m-1,4);
    
    % let's try Collard's initialization
    v0 = (((xgrid.^alph).^(1-sig)-1)/((1-sig)*(1-bet))); % you get exactly the same thing
    
    % matlab's spline
    pp = spline(xgrid,v0);
    fake_kp = (xgrid(2)-xgrid(1))/2;
    vhat_m = ppval(pp,fake_kp);
    
    % return
    
    % Let's try out the spline functions
    % v0_vector = vhat_spline(coeffs0, xgrid,xgrid);
    % v0_scalar = vhat_spline(coeffs0,xgrid, xgrid(end-1));
    %
    % v_try_scalar = -TV_spline(xgrid(1),coeffs0,xgrid(1),xgrid, param);
    
    % Let's try my spline on simple things
    % x^2
    % interval = linspace(-1,1,20);
    % % Approximand in "dataset"
    % true_fun = @(x) x.^2;
    % f = true_fun(interval);
    %
    % objh = @(coeffs) obj_spline_secant_hermite(coeffs,interval,f);
    %     [coeffs1,FVAL,exitflag] = fsolve(objh,coeffs0, options2);
    % fhat = vhat_spline(coeffs1,interval,interval);
    % fhat_end = vhat_spline(coeffs1,interval,interval(end));
    % figure
    % plot(f); hold on
    % plot(fhat)
    % yo it seems to work so nicely
    
    % % four known points
    % xi = [-1,0,1,2];
    % yi = [2,0,-2,0]';
    % M = [1,-1,1,-1; 1,0,0,0; 1,1,1,1; 1,2,4,8];
    % a = M^(-1)*yi;
    % a0 = zeros(size(a));
    % objh = @(a) obj_spline_secant_hermite(a,xi,yi);
    % [a1,FVAL,exitflag(iter)] = fsolve(objh,a0', options2); % my approach doesn't work here at all
    
    % return
    
    % coeffs=coeffs0;
    coeffs=pp;
    v=v0;
    
    v1 = zeros(m,1);
    kp = zeros(m,1);
    exitflag=nan;
    % always start searching at the beginning of the grid
    k0 = xgrid(1);
    datestr(now)
    tic
    while crit > epsi && iter< maxiter
        % Step 1: maximization
        for j=1:m
            kj=xgrid(j);
            negv = @(kp) TV_spline(kp,coeffs,kj,xgrid, param);
            [kp(j,iter)] = fminunc(negv, k0, options1);
            
            % compute the value function at the maximizing kp
            v1(j) = -TV_spline(kp(j,iter),coeffs,kj, xgrid, param);
        end
        
        % Step 2: fitting
        %     objh = @(coeffs) obj_spline_secant_hermite(coeffs,xgrid,v1);
        %     % I think this is the conceptually correct thing: the new k are the nodes; or you know, maybe not, maybe the nodes for the spline should be the same %
        %     %     objh = @(coeffs) obj_spline_secant_hermite(coeffs,kp(:,iter),v1);
        %     [coeffs1,FVAL,exitflag(iter)] = fsolve(objh,coeffs0, options2);
        %     %     %This is awful: it's totally diverging
        %     %     objh = @(coeffs)obj_vhat_spline(coeffs,xgrid,v1);
        %     %     [coeffs1,FVAL,exitflag(iter)] = fsolve(objh,coeffs, options2);
        
        % this works... juhuuuu
        coeffs1 = spline(xgrid,v1);
        %     return
        % Compute stopping criterion and update
        crit = max(abs(v1-v));
        if mod(iter,50)==0
            disp(['Concluding iteration = ', num2str(iter)])
            disp(['Stopping criterion = ', num2str(crit)])
        end
        iter=iter+1;
        coeffs=coeffs1;
        v=v1;
    end
    toc
    
    % Compute consumption function
    c = xgrid.^alph+(1-delt).*xgrid-kp(:,end);
    
    
    figure
    plot(xgrid,v)
    title('Collard''s Fig. 7.7 Value function against k_t - SPLINE')
    figname = 'vfi_spline_valuefun';
    
    if print_figs ==1
        disp(figname)
        cd(figpath)
        export_fig(figname)
        cd(current_dir)
        close
    end
    
    
    figure
    subplot(1,2,1)
    plot(xgrid,kp(:,end)); hold on
    plot(xgrid,xgrid,'k--')
    title('k_{t+1} as a function of k_t')
    subplot(1,2,2)
    plot(xgrid,c)
    title('c_{t} as a function of k_t')
    sgtitle('Collard''s Fig. 7.8 Decision rules - SPLINE')
    figname = 'vfi_spline_decisions';
    
    if print_figs ==1
        disp(figname)
        cd(figpath)
        export_fig(figname)
        cd(current_dir)
        close
    end
    
end

%% let's try parametric VFI with a Piecewise Cubic Hermite Interpolating Polynomial. PCHIP

clc

sig  = 1.5;
delt = 0.1;
bet  = 0.95;
alph = 0.3;
param = [sig,delt,bet,alph];

%Optimization Parameters
% options = optimset('fmincon');
options1 = optimset('fminunc');
options1 = optimset(options1, 'TolFun', 1e-16, 'display', 'none');

m=20; % number of nodes (index i), m-1 intervals
crit=1;
iter=1;
maxiter=600;%300
epsi=1e-6; %-6
ks = ((1-bet*(1-delt))/(alph*bet))^(1/(alph-1));
dev = 0.9;
kmin = (1-dev)*ks;
kmax = (1+dev)*ks;
% for the pchip, uniformly placed nodes are fine
xgrid = linspace(kmin,kmax,m)';

% let's try Collard's initialization
v0 = (((xgrid.^alph).^(1-sig)-1)/((1-sig)*(1-bet))); % you get exactly the same thing

% matlab's spline
pp = pchip(xgrid,v0);
fake_kp = (xgrid(2)-xgrid(1))/2;
vhat_m = ppval(pp,fake_kp);

% return

coeffs=pp;
v=v0;

v1 = zeros(m,1);
kp = zeros(m,1);
exitflag=nan;
% always start searching at the beginning of the grid
k0 = xgrid(1);
datestr(now)
tic
while crit > epsi && iter< maxiter
    % Step 1: maximization
    for j=1:m
        kj=xgrid(j);
        negv = @(kp) TV_pchip(kp,coeffs,kj, param);
        [kp(j,iter)] = fminunc(negv, k0, options1);
        
        % compute the value function at the maximizing kp
        v1(j) = -TV_pchip(kp(j,iter),coeffs,kj, param);
    end
    
    % Step 2: fitting
    coeffs1 = pchip(xgrid,v1);
    
    % Compute stopping criterion and update
    crit = max(abs(v1-v));
    if mod(iter,50)==0
        disp(['Concluding iteration = ', num2str(iter)])
        disp(['Stopping criterion = ', num2str(crit)])
    end
    iter=iter+1;
    coeffs=coeffs1;
    v=v1;
end
toc

% Compute consumption function
c = xgrid.^alph+(1-delt).*xgrid-kp(:,end);


figure
plot(xgrid,v)
title('Collard''s Fig. 7.7 Value function against k_t - PCHIP')
figname = 'vfi_pchip_valuefun';

if print_figs ==1
    disp(figname)
    cd(figpath)
    export_fig(figname)
    cd(current_dir)
    close
end


figure
subplot(1,2,1)
plot(xgrid,kp(:,end)); hold on
plot(xgrid,xgrid,'k--')
title('k_{t+1} as a function of k_t')
subplot(1,2,2)
plot(xgrid,c)
title('c_{t} as a function of k_t')
sgtitle('Collard''s Fig. 7.8 Decision rules - PCHIP')
figname = 'vfi_pchip_decisions';

if print_figs ==1
    disp(figname)
    cd(figpath)
    export_fig(figname)
    cd(current_dir)
    close
end

%% Stochastic value function iteration - Collard's code

sigma   = 1.50;
delta   = 0.10;
beta    = 0.95;
alpha   = 0.30;
rho     = 0.80;
se      = 0.12; % volatility of the shock
nbk = 1000; % number of kapital gridpoints
nba =2;  % number of technology gridpoints
crit =1; 
epsi = 1e-6; 
% Discretization of the shock - Collard does this w/o Tauchen-Hussey b/c
% tech shock only takes on 2 values
p       = (1+rho)/2;
PI      = [p 1-p;1-p p];
se = 0.12; 
ab =0;
a1      = exp(-se*se/(1-rho*rho));
a2      = exp(se*se/(1-rho*rho));
A = [a1 a2];
% Discretization of the state space - just a grid
kmin = 0.2;
kmax =6;
kgrid = linspace(kmin,kmax,nbk)'; 
c = zeros(nbk,nba);
util = zeros(nbk,nba);
v = zeros(nbk,nba);
Tv = zeros(nbk,nba);
% Collard doesn't seem to define k anywhere, I think it's just the kgrid
k= kgrid;
while crit>epsi;
    for i=1:nbk
        for j=1:nba;
            c           = A(j)*k(i)^alpha+(1-delta)*k(i)-k;
            neg         = find(c<=0);
            c(neg)      = NaN;
            util(:,j)   = (c.^(1-sigma)-1)/(1-sigma);
            util(neg,j) = -1e12;
        end
        [Tv(i,:),dr(i,:)] = max(util+beta*(v*PI));
    end;
    crit = max(max(abs(Tv-v)));
    v        = Tv;
    iter = iter+1;
end

kp      = k(dr);
for j=1:nba;
    c(:,j)   = A(j)*k.^alpha+(1-delta)*k-kp(:,j);
end

figure
plot(k,v(:,1)); hold on
plot(k,v(:,2))
title('Collard''s Fig. 7.9 Value function against k_t - with bivariate tech shock, discretized')
figname = 'vfi_stoch_valuefun';

if print_figs ==1
    disp(figname)
    cd(figpath)
    export_fig(figname)
    cd(current_dir)
    close
end


figure
subplot(1,2,1)
plot(k,kp(:,1)); hold on
plot(k,kp(:,2))
plot(k,k,'k--')
title('k_{t+1} as a function of k_t')
subplot(1,2,2)
plot(k,c(:,1)); hold on
plot(k,c(:,2))
title('c_{t} as a function of k_t')
sgtitle('Collard''s Fig. 7.10 Decision rules - with bivariate tech shock, discretized')
figname = 'vfi_stoch_decisions';

if print_figs ==1
    disp(figname)
    cd(figpath)
    export_fig(figname)
    cd(current_dir)
    close
end

%% Let me do my own functions

function Tx = chebyshev(x,n) % checked - correct
[xm,xn]=size(x);
tau=max(xm,xn); % the number of data points
x=reshape(x,tau,1);
Tx = zeros(tau,n+1);
% creates Chebyshev polynomials for all x, of order 0,...,n
Tx(:,1) = ones(tau,1); % order 0
Tx(:,2) = x.*ones(tau,1);
for j=3:n+1
    Tx(:,j) = 2*x.*Tx(:,j-1)-Tx(:,j-2);
end
end

function vhat = vhat_cheby(b,x,n,xmin,xmax) % checked - correct
% Make sure b has the right size
b = reshape(b,n+1,1);
% Convert x back to z
z = 2*(x-xmin)/(xmax-xmin)-1;
% Evaluate Chebyshevs at z
T = chebyshev(z,n);
% Evaluate the approximated value function given the coefficients
vhat = T*b;
end

function negvj = TV(kp,bi,xj,n,xmin,xmax, param)
% this function just gives you today's value for this state xj, times -1 (it's not a min)
sig = param(1); delt=param(2); bet=param(3); alph=param(4);
% we're moving along different values for today's capital, and we wanna
% choose tomorrow's optimally, kp
k=xj;
% compute c given k, util given c
kp = sqrt(kp.^2); % insures positivity of k' - is this necessary? I suspect not.
c = k.^alph+(1-delt)*k-kp; % computes consumption
d = find(c<=0);% find negative consumption
c(d) = NaN; % get rid of negative c
util = (c.^(1-sig)-1)/(1-sig); % computes utility
util(d) = -1e12; % utility = low number for c<0
% compute value function given k
v = vhat_cheby(bi,kp,n,xmin,xmax);
% finally, the T-map TVhat:
vj = util + bet*v;
negvj =-vj;
end

function negvj = TV_spline(kp,coeffs,k,xgrid, param)
% this function just gives you today's value for this state xj, times -1 (it's not a min)
sig = param(1); delt=param(2); bet=param(3); alph=param(4);
% we're moving along different values for today's capital, and we wanna
% choose tomorrow's optimally, kp
% compute c given k, util given c
kp = sqrt(kp.^2); % insures positivity of k' - is this necessary? I suspect not.
c = k.^alph+(1-delt)*k-kp; % computes consumption
d = find(c<=0);% find negative consumption
c(d) = NaN; % get rid of negative c
util = (c.^(1-sig)-1)/(1-sig); % computes utility
util(d) = -1e12; % utility = low number for c<0
% compute value function given k
% v = vhat_spline(coeffs,xgrid,kp);
% try matlab's spline
v = ppval(coeffs,kp);

% finally, the T-map TVhat:
vj = util + bet*v;
negvj =-vj;
end


function negvj = TV_pchip(kp,coeffs,k, param)
% this function just gives you today's value for this state xj, times -1 (it's not a min)
sig = param(1); delt=param(2); bet=param(3); alph=param(4);
% we're moving along different values for today's capital, and we wanna
% choose tomorrow's optimally, kp
% compute c given k, util given c
kp = sqrt(kp.^2); % insures positivity of k' - is this necessary? I suspect not.
c = k.^alph+(1-delt)*k-kp; % computes consumption
d = find(c<=0);% find negative consumption
c(d) = NaN; % get rid of negative c
util = (c.^(1-sig)-1)/(1-sig); % computes utility
util(d) = -1e12; % utility = low number for c<0
% compute value function given k
v = ppval(coeffs,kp);

% finally, the T-map TVhat:
vj = util + bet*v;
negvj =-vj;
end



%%%%%%%%%%%%%%%%%%%%%%%  RETIRED %%%%%%%%%%%%%%%%%%%%%%%
function vhat = vhat_spline(coeffs, xgrid, x)
ngrid=length(xgrid);
T = length(x);
% fish out coefficients
a = coeffs(:,1);
b = coeffs(:,2);
c = coeffs(:,3);
d = coeffs(:,4);

sx = zeros(T,1);
for i=1:T
    % find the first node larger than x
    next_node = find(xgrid>x(i),1);
    int_i = next_node-1;
    if int_i==0
        int_i=1;
    elseif x(i)>=max(xgrid)
        int_i=ngrid-1;
    end
    sx(i) = a(int_i)+ b(int_i)*x(i) + c(int_i)*x(i)^2 + d(int_i)*x(i)^3;
end
vhat = sx;
end

function resids = obj_vhat_spline(a,xgrid,v)
vhat = vhat_spline(a,xgrid,xgrid);
resids = (vhat-v).^2;
end

function b = compute_cheby_coeffs(v,x,n,xmin,xmax) % not used b/c not necessary!
b=zeros(n+1,1);
% convert xi to zi
z = (2.*x -xmin -xmax)./(xmax-xmin);
[zm,zn] = size(z);
m = max(zm,zn); % recover the number of nodes
% make sure size of v is right
v=reshape(v,1,m);
% compute chebyshevs
T = chebyshev(z,n);
b(1) = 1/m*(sum(v));
for j=2:n+1
    viTj = v*T(j,:)';
    b(j) = 2/m * viTj;
end

end
%% Collard's functions
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%% A nice thing is Collard's intuitive Chebyshev and value functions
% chebyshev = @(Tn,Tn_1,n,x) 2.*x.*Tn -Tn_1; % this was mine and it's not as convenient as Collard's
% Tx =chebychev([-1:0.01:1],5); % ok it's working
% plot(Tx(:,5))

function Tx=chebychev(x,n)
X=x(:);
lx=size(X,1);
if n<0
    error('n should be a positive integer');
end
switch n
    case 0
        Tx=[ones(lx,1)];
    case 1
        Tx=[ones(lx,1) X];
    otherwise
        Tx=[ones(lx,1) X];
        for i=3:n+1
            Tx=[Tx 2*X.*Tx(:,i-1)-Tx(:,i-2)];
        end
end
end

% the other nice thing is you can see how the value function is a function
% given parameters
function v = value(k,param,theta)
kmin    = param(1);
kmax    = param(2);
n       = param(3);
zk      = 2*(k-kmin)/(kmax-kmin)-1; % convert k back to zk of [-1,1]
v       = chebychev(zk,n)*theta; % this is the polynomial approx: chebyshev polynomial * coefficients
end

% and finally the T-map, mapping V(k') to V
function res=tv(kp,param,theta)
alpha   = param(1);
beta    = param(2);
delta   = param(3);
sigma   = param(4);
kmin    = param(5);
kmax    = param(6);
n       = param(7);
k       = param(8);
kp      = sqrt(kp.^2); % insures positivity of k?
v       = value(kp,[kmin kmax n],theta); % computes the value function
c       = k.^alpha+(1-delta)*k-kp; % computes consumption
d       = find(c<=0);% find negative consumption
c(d)    = NaN; % get rid of negative c
util    = (c.^(1-sigma)-1)/(1-sigma); % computes utility
util(d) = -1e12; % utility = low number for c<0
res     = -(util+beta*v); % compute -TV (we are minimizing)
end