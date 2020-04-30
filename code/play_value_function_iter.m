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
sig  = 1.5;
delt = 0.1;
bet  = 0.95;
alph = 0.3;

ngrid=20; % Collard's nbk
p = 10; % order of Chebyshev polyonomials (0,1,...,10)
crit=1;
iter=1;
epsi=1e-6;
ks = ((1-bet*(1-delt))/(alph*bet))^(1/(alph-1));
dev = 0.9;
kmin = (1-dev)*ks;
kmax = (1+dev)*ks;
zk = -cos(2*[1:ngrid]'-1)*pi/(2*(ngrid));
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
while crit>epsi
    k0 = kgrid(1);
    for i=1:ngrid
        param = [alph, bet, delt, sig, kmin, kmax, ngrid, kgrid(i)];
        % kp(i) = some optimization routine that gives you optimal k'
        % Given current coefficients th0 and today's k, kgrid(i), starting
        % with a guess of k0, solve for k':
        objh = @(kp) tv(kp,param,th0);
        [kp(i)] = fmincon(objh, k0, [],[],[],[],kmin,kmax,[],options);
        % fmincon(FUN,X0,A,B,Aeq,Beq,LB,UB,NONLCON,OPTIONS)
        % I wonder if we could get around doing fmincon...
        
        k0 = kp(i);
        %get back today's value function given the optimal kp(i) and the
        %T-map given current "estimated" coefficients theta0:
        Tv(i) = -tv(kp(i), param, th0);
    end
    % Fitting step: Given value function today, back out coefficients
    thet =X\Tv;
    crit = max(abs(Tv-v));
    v = Tv;
    th0 = thet;
    iter=iter+1;
end
toc % indeed it takes 242 iterations (but much longer than the above: 60 seconds)

% Compute consumption function
kgrid
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
v       = chebychev(zk,n)*theta; % this is the polyonomial approx: chebyshev polynomial * coefficients
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
