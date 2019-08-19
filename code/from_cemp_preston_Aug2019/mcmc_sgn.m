function [param1,p1,accept] = mcmc_sgn(param0,Y,postold,step,lbnd,ubnd,pshape,pmean,pstdd,N, obs_blocksize)
% Is actually the MCMC algorithm for CEMP for the case of no missings in
% the data.
% pshape = what PDF does the prior have (this and the next two are (nparam x 1)
% pmean = prior mean
% pstdd = prior std dev.

accept = 0;
nparam = length(param0);
T = size(Y,2);

% Adapt VC-matrix of proposal
if size(step,1) == 1 % if step is diag(VC_hat) (a vector)
    VC_p = step.*eye(nparam); % VC-matrix of the proposal.
elseif size(step,1) >= 2 % if step IS VC_hat
    VC_p = step;
else
    warning('Stepsize incorrect.')
end

% %Step 1: Draw a new set of parameters from proposal (here uniform) distribution
% fstep = 2*(rand(1,nparam)-.5); %Point from -1 to 1 of step size. Draw from Unif(-1,1).
% paramp = param0+fstep.*step;   %Theta': the new draw of parameters. Center them around theta^t, and scale the draw by stepsize step.

% Do step 1, drawing from a normal
% paramp = param0 + step.*randn(1,nparam);
paramp = param0 + (VC_p*randn(nparam,1))';
% paramp = mvnrnd(param0,VC_p,1);

% % Do step 1, drawing from a Student's t % doesn't work yet
% V = T-1; % degrees of freedom
% t = trnd(V, 1,nparam);

%Step 1.2: Evaluate prior distribution for proposed paramp (this draws on Dongho, RandomWalkMH.m and objfcn.m)
lnprior = 0; % this is actually the prior probability and it's obtained by adding all the prior probs for all params up.
for i = 1:nparam
    if pshape(i) == 1      % BETA Prior
        a = (1-pmean(i))*pmean(i)^2/pstdd(i)^2 - pmean(i);
        b = a*(1/pmean(i) - 1);
        lnprior = lnprior + log( pdf('beta',paramp(i),a,b) );
        
    elseif pshape(i) == 2  % GAMMA PRIOR
        b = pstdd(i)^2/pmean(i);
        a = pmean(i)/b;
        lnprior = lnprior + log( pdf('gam',paramp(i),a,b) );
        
    elseif pshape(i) == 3  % GAUSSIAN PRIOR
        a = pmean(i);
        b = pstdd(i);
        lnprior = lnprior + log( pdf('norm',paramp(i),a,b) );
        
    elseif pshape(i) == 4 % INV WISHART PRIOR
        % this is not mean nor stdev
        a = pmean(i);
        b = pstdd(i);
        lnprior = lnprior + (b/2)*log(det(a)) - ((b+2)/2)*log(det(paramp(i))) - 0.5*trace(a*inv(paramp(i)));
        
    elseif pshape(i) == 5  % uniform prior
        a = pmean(i);
        b = pstdd(i);
        lnprior = lnprior + log(1/(b-a)*all(paramp(i)>a & paramp(i)<b));
    elseif pshape(i) == 6  % INV GAMMA PRIOR
        aux = (pstdd(i)/pmean(i))^2;
        a = (1-2*aux)/(1-aux);
        b = pmean(i)*(a-1);
        lnprior = lnprior + log( b^a/gamma(a)*paramp(i).^(-a-1).*exp(-b/paramp(i)) );
    end
end


%Step 2: Compute liklihood at new set of parameters, and compare
if sum(paramp<lbnd) + sum(paramp>ubnd)
    b = 0;  %Reject if outside of bnds (aka priors)
    param1 = param0;
    p1 = postold;
    return;
else
    %%% COMPARE MPFs here
    %     pp = MPF_test_resampling1(Y,paramp,N,obs_blocksize); % <--- default (16 July 2019)
    % % %     pp = MPF_test_resampling2(Y,paramp,N,obs_blocksize);
    % estimate params one by one
    %     pp = MPF_test_resampling1_1by1(Y,paramp,N,obs_blocksize);
    gbar = param0(3);
    randgbar = gbar.*rand(1,N); % CEMP keep this constant
    k10 = randgbar.^(-1); % CEMP keep this constant
    pibar10 = randn(1,N); % CEMP keep this constant
    predrand = randn(T,N); % generate all the random draws for the nonlinear prediction step outside the filter
    pp = MPFquick(Y,param0,N,k10,pibar10,predrand,obs_blocksize);
    
    
    % I added the below line: posterior is log(likelihood) + log(prior)
    postnew = pp + lnprior;
    b = exp(postnew-postold);  %Ratio of liklihoods
end


%Step 3: Choose to accept or reject
if b >= 1 %Accept for sure
    param1 = paramp;
    p1 = pp;
    accept = 1;
elseif rand(1,1)<=b %Accept w/p b;
    param1 = paramp;
    p1 = pp;
    accept=1;
else
    param1 = param0;
    p1 = postold;
end