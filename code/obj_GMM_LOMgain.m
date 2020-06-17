function loss = obj_GMM_LOMgain(alph,x,xxgrid, yygrid,param,gx,hx,eta,e,T,ndrop,PLM,gain,Om_data, W1)
% alph are the coefficients, x is the grid
% 9 June 2020

% check "global" nonnegativity of k1
k10 = ndim_simplex_eval(x,[xxgrid(:)';yygrid(:)'],alph);
if min(k10)<0
    loss = 1e+10;
    disp('k1 was negative on fine grid, not even bothering to do simulation')
else
    
    [ny, nx] = size(gx);
    
    knowTR=1;
    mpshock=1;
    % Simulate data given parameters
    [~, y, k] = sim_learnLH_clean_approx(alph,x,param,gx,hx,eta, PLM, gain, T+ndrop,ndrop,e, knowTR,mpshock);
    
    
    k1 = 1./k(1:end-1); % cut off last period where k is unset
    % Do not filter data and estimate VARs if the current coefficients
    % alpha lead to an explosive learning simulation
    if isinf(max(k1)) || min(k1)<0
        loss = 1e+10;
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
        ystar = nan(ny,T-2*K);
        for i=1:ny
            ystar(i,:) = BKfilter(y(i,:)');
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
        max_lags = 6;
        % [AIC,BIC,HQ] = aic_bic_hq(filt',max_lags);
        % p =min([AIC,BIC,HQ]);
        p = 1;
        % A is the impact matrix, identified via Cholesky, B is the beta_ols, res are
        % the residuals, sigma is the estimated VC matrix.
        [A,B,res,sigma] = sr_var(filt', p);
        
        % Rewrite the VAR(p) as VAR(1) (see Hamilton, p. 273, Mac)
        PHI = B(2:end,:)';
        F = [PHI; [eye(ny*(p-1)), zeros(ny*(p-1),ny)]];
        Q = [[sigma, zeros(ny,ny*(p-1))]; zeros(ny*(p-1),ny*p)];
        % check sizes
        np = ny*p;
        
        vecSig = (eye(np^2)-kron(F,F))\vec(Q);
        % VC matrix of data y
        Gamj = zeros(ny,ny,K+1);
        Sig = reshape(vecSig,np,np);
        for j=0:K
            % jth Autocov of data y, still as a VAR(1)
            Sigj = F^j * Sig;
            % jth Autocov of data y, finally as a VAR(p)
            Gamj(:,:,j+1) = Sigj(1:ny,1:ny);
        end
        % moments vector
        Om = vec(Gamj);
        
        % Compute GMM loss
        loss = (Om_data -Om)'*W1*(Om_data - Om);
    end
end