function lnprior = evaluate_prior_materials3(param,pshape,pmean,pstdd)
nparam = length(param);
lnprior = 0; % this is actually the prior probability and it's obtained by adding all the prior probs for all params up. (1x1)
for i = 1:nparam
    if pshape(i) == 1      % BETA Prior
        a = (1-pmean(i))*pmean(i)^2/pstdd(i)^2 - pmean(i);
        b = a*(1/pmean(i) - 1);
        lnprior = lnprior + log( pdf('beta',param(i),a,b) );
        
    elseif pshape(i) == 2  % GAMMA PRIOR
        b = pstdd(i)^2/pmean(i);
        a = pmean(i)/b;
        lnprior = lnprior + log( pdf('gam',param(i),a,b) );
        
    elseif pshape(i) == 3  % GAUSSIAN PRIOR
        a = pmean(i);
        b = pstdd(i);
        lnprior = lnprior + log( pdf('norm',param(i),a,b) );
        
    elseif pshape(i) == 4 % INV WISHART PRIOR
        % this is not mean nor stdev
        a = pmean(i);
        b = pstdd(i);
        lnprior = lnprior + (b/2)*log(det(a)) - ((b+2)/2)*log(det(param(i))) - 0.5*trace(a*inv(param(i)));
        
    elseif pshape(i) == 5  % uniform prior
        a = pmean(i);
        b = pstdd(i);
        lnprior = lnprior + log(1/(b-a)*all(param(i)>a & param(i)<b));
    elseif pshape(i) == 6  % INV GAMMA PRIOR
        aux = (pstdd(i)/pmean(i))^2;
        a = (1-2*aux)/(1-aux);
        b = pmean(i)*(a-1);
        lnprior = lnprior + log( b^a/gamma(a)*param(i).^(-a-1).*exp(-b/param(i)) );
    end
end