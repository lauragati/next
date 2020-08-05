function [res, Om, FEt_1] = obj_GMM_LOMgain_univariate(alph,x,xxgrid,param,gx,hx,eta,e,v,T,ndrop,PLM,gain,p,Om_data, W1,Wdiffs2,Wmid,Wmean,use_expectations_data,alph0,Wprior)
% alph are the coefficients, x is the grid
% 9 June 2020
% Update 17 June 2020: rewritten to work with lsqnonlin
% Update 24 June 2020: no longer need to input alph0 and Wprior if you're
% not penalizing deviations from the prior.
this_code = mfilename;
max_no_inputs = nargin(this_code);
if nargin < max_no_inputs
    Wprior = 0;
    alph0 = 0;
end

% disp('Current guess alpha = ')
% disp(num2str(alph))
% check "global" nonnegativity of k1
k10 = ndim_simplex_eval(x,xxgrid(:)',alph);
if min(k10)<0
    res = 1e+10*ones(size(Om_data));
    disp('k1 was negative on fine grid, not even bothering to do simulation')
else
        
    knowTR=1;
    mpshock=1;
    
    % Simulate data given parameters
    [~, y, k,phi,~,~,FEt_1] = sim_learnLH_clean_approx_univariate(alph,x,param,gx,hx,eta, PLM, gain, T+ndrop,ndrop,e,v, knowTR,mpshock);
    
    if use_expectations_data == 0
        y_data = y(:,1:end-1);
    else     % add one-step forecast to data
        y_data = [y(:,1:end-1); squeeze(phi(1,1,1:end-1))'];
    end
    nobs = size(y_data,1);

    k1 = 1./k(1:end-1); % cut off last period where k is unset
    % Do not filter data and estimate VARs if the current coefficients
    % alpha lead to an explosive learning simulation
    if isinf(max(k1)) || min(k1)<0
        res = 1e+10*ones(size(Om_data));
        Om =nan;
    else
        % Filter the simulated data
        % % 1) HP filter
        % g = nan(size(y));
        % c = nan(size(y));
        % for i=1:ny
        %     [g(i,:),c(i,:)] = HPfilter(y(i,:)');
        % end
        
        % 2) Hamilton filter
        % h=8;
        % v = nan(ny,T-4-h+1);
        % for i=1:ny
        %     [v(i,:)] = Hamiltonfilter(y(i,:)');
        % end
        
        % 3) BK filter
        K=12;
        ystar = nan(nobs,T-2*K-1);
        for i=1:nobs
            ystar(i,:) = BKfilter(y_data(i,:)');
        end
        
        filt=ystar;
        ranky = rank(ystar'*ystar);
        %     if ranky < ny
        %         warning('Model-generated data matrix is not full rank')
        %     end
        
        % Compute the model-implied moments
        % compute moments, Om, as the autocovariances of the data for lags
        % 0,1,...,K
        K=4;
        % Take the initial data, estimate a VAR
        % using the same lags p, K as for the real data
        [~,B,~,sigma] = sr_var(filt', p);
        
        % Rewrite the VAR(p) as VAR(1) (see Hamilton, p. 273, Mac)
        PHI = B(2:end,:)';
        F = [PHI; [eye(nobs*(p-1)), zeros(nobs*(p-1),nobs)]];
        Q = [[sigma, zeros(nobs,nobs*(p-1))]; zeros(nobs*(p-1),nobs*p)];
        % check sizes
        np = nobs*p;
        
        vecSig = (eye(np^2)-kron(F,F))\vec(Q); % Hamilton p. 279 Mac, eq. [10.2.18]
        % VC matrix of data y
        Gamj = zeros(nobs,nobs,K+1);
        Gamj_own = zeros(nobs,K+1);
        
        Sig = reshape(vecSig,np,np);
        for j=0:K
            % jth Autocov of data y, still as a VAR(1)
            Sigj = F^j * Sig;
            % jth Autocov of data y, finally as a VAR(p)
            Gamj(:,:,j+1) = Sigj(1:nobs,1:nobs);
            % This only takes the covariances wrt the same var: (e.g Cov(x_t,x_{t-1}))
            Gamj_own(:,j+1) = diag(Sigj(1:nobs,1:nobs));
        end
        % moments vector
        Om = vec(Gamj);
        %         Om = vec(Gamj_own);
        
        % additional moments
        % 15 July 2020 correction for convexity moment (see Notes)
        % First derivatives:
        deltalph  = diff(alph,1);
        fegrid = x{1};
        deltfegrid = diff(fegrid,1);
        yp = deltalph./deltfegrid';
        
        % Second derivatives:
        deltyp =  diff(yp,1);
        ypp    = deltyp ./ deltfegrid(2:end)';
        
        nfx2 = length(ypp);
        diffs2_moment = (ypp.*(ypp<=0)).^2;
        calibrated_moment = (mean(k1)-0.05)^2;
        nfe = length(alph);
        if mod(nfe,2)==0 % if there are even no. of knots, nevermind
            alph_mid = 0;
        else
            midpoint=ceil(nfe/2);
            alph_mid = alph(midpoint);
        end
        
        
        % Compute GMM loss, not squared, just weighted ("weighted, not squared")
        devprior = sum(abs(alph0-alph))*Wprior;
        res = (Om_data -Om).*diag(W1);
        % add extra moments
        res(end+1) = devprior;
        res(end+1) = Wmean*calibrated_moment;
        res(end+1) = Wmid*alph_mid;
        res(end+1:end+nfx2) = Wdiffs2*diffs2_moment;
        
    end
end