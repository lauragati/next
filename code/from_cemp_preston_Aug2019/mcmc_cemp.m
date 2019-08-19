function [param1,p1,accept] = mcmc_cemp(param0,nanY,Y1,Y2,Y3,postold,step,lbnd,ubnd,pshape,pmean,pstdd,N)
% 1 July 2019
% pmean = prior mean
% pstdd = prior std dev.

accept = 0;
nparam = length(param0);

% %Step 1: Draw a new set of parameters from proposal (here uniform) distribution
% fstep = 2*(rand(1,nparam)-.5); %Point from -1 to 1 of step size. Draw from Unif(-1,1).
% paramp = param0+fstep.*step;   %Theta': the new draw of parameters. Center them around theta^t, and scale the draw by stepsize step.

% Do step 1, drawing from a normal
paramp = param0 + step.*randn(1,nparam);

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
    pp = MPF_test_missings(nanY,Y1,Y2,Y3,paramp,N);

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