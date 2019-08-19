function LL = particlefilter_generic(Y,k,m)
% Uses notation of Kitagawa (1996), p. 6, top.
% See materials3.m of experimentation for the original implementation.
% There I've also updated the resampling, which ISN'T updated here. 
% Y = data: (l,T)
% k = number of states
% m = number of particles

% The system looks like:
% xn = F(x_{n-1}, vn)
% yn = H(xn, wn)
% where n indexes time, x is (k,1), y is (l,1)

% You'll need to copy this general structure and modify it every time for specific models.
% In particular, specify model-specific functions F_fun, G_fun and constants dGdyn.

% Also, change around distributions/densities as needed.
% -------------------------------------------------------
T = size(Y,2);
% 1.) Generate f0^j ~ p0(x), j=1,...,m.
f0 = randn(k,m); % I'm assuming I know somehow that p0(x) = N(0,1).
fn_1 =f0; % initialization
% 2.)
p = zeros(k,m); % these are the particles
alph = zeros(T,m); % likelihoods of the particles
fn = zeros(size(fn_1));

for n=1:T
    % a) Generate vn^j ~ q(v) j=1,...,m. 
    vn = randn(k,m);
    % b) Generate pn^j = F(f(n-1)^j, vn^j) for j=1,...,m.
    for j=1:m
        p(:,j) =  F_fun(fn_1(:,j),vn(:,m));
    % c) Compute alpha_n^j = r( G(yn,pn^j) ) * |dG/dyn| for j=1,...,m.
        G = G_fun(Y(:,n), p(:,j));
        MU = 4; % <---- specify this
        SIG = eye(k); % <---- specify this (VC matrix of measurement error)
        r = @(x) mvnpdf(x,MU,SIG); % density of observation noise w. 
        alph(n,j) = r(G)* dGdyn;
    end
    % d) Resampling: fn^j for j=1,...,m.
    prob_weights = zeros(m,1);
    for j=1:m
        if isnan(alph(n,j)/sum(alph(n,:)))
            prob_weights(j) = 0;
        else
            prob_weights(j) = alph(n,j)/sum(alph(n,:));
        end
    end
    for j=1:m
        index = randsample(1:m,1,true,prob_weights);
        fn(:,j) = p(n,:,index);
    end
    fn_1 = fn; % update simulated states
end

LL = -T*log(m);
for n=1:T
    LL = LL + sum(alph(n,:));
end