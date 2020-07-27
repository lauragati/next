function [res, Om, Om_n] = obj_GMM_check_mean(alph,x,rho,eN,T,ndrop,p,Om_data,W1)

[~,~,N] = size(eN);

% states and jumps
Om_n = nan(numel(Om_data),N);

for n=1:N
    s = zeros(T+ndrop,1);
    y = zeros(T+ndrop,1);

    e=squeeze(eN(:,:,n));
    s(1,n) = e(1);
    y(1,n) = 0;
    for t=2:T+ndrop
        s(t,n)=rho*s(t-1) +e(t);
        y(t,n)= ndim_simplex_eval(x, s(t),alph);
    end
    
    % Map into command_acf_sim_data_univariate.m
    
    y_data = y(ndrop+1:end)';
    nobs = size(y_data,1);
    % Filter the simulated data
    % 1) HP filter
    g = nan(size(y_data));
    c = nan(size(y_data));
    for i=1:nobs
        [g(i,:),c(i,:)] = HPfilter(y_data(i,:)');
    end
    
    
    % 2) Hamilton filter
    % h=8;
    % v = nan(ny,T-4-h+1);
    % for i=1:ny
    %     [v(i,:)] = Hamiltonfilter(y(i,:)');
    % end
    
%     % 3) BK filter
%     K=12;
%     ystar = nan(nobs,T-2*K);
%     for i=1:nobs
%         ystar(i,:) = BKfilter(y_data(i,:)');
%     end
    
    filt=c;
    ranky = rank(filt'*filt);
        if ranky < nobs
            warning('Model-generated data matrix is not full rank')
        end
    
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
    
    vecSig = (eye(np^2)-kron(F,F))\vec(Q);
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
    Om_n(:,n) = vec(Gamj);
    %         Om = vec(Gamj_own);
end
Om = nanmean(Om_n,2);

res = (Om_data -Om).*diag(W1);


end