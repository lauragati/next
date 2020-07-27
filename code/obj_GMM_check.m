function [res, Om] = obj_GMM_check(alph,x,rho,e,T,ndrop,p,Om_data,W1);
% states and jumps
s = zeros(T+ndrop,1);
y = zeros(T+ndrop,1);

s(1) = e(1);
y(1) = 0;
for t=2:T+ndrop
    s(t)=rho*s(t-1) +e(t);
    y(t)= ndim_simplex_eval(x, s(t),alph);
end

y = y(ndrop+1:end);

% Map into command_acf_sim_data_univariate.m
y_data = y';
[nobs,~] = size(y_data);

% Filter the data
% HP filter
g_data = nan(size(y_data));
c_data = nan(size(y_data));
for i=1:nobs
    [g_data(i,:),c_data(i,:)] = HPfilter(y_data(i,:)');
end
% 
% % Hamilton filter
% h=8;
% v_data = nan(nobs,T-4-h+1);
% for i=1:nobs
%     [v_data(i,:)] = Hamiltonfilter(y_data(i,:)');
% end
% lost_periods_H = h+3;


% % BK filter
% K=12;
% ystar_data = nan(nobs,T-2*K);
% for i=1:nobs
%     ystar_data(i,:) = BKfilter(y_data(i,:)');
% end
% lost_periods_BK = 2*K;


% choose your favorite filtered data
filt_data = c_data;
% lost_periods = lost_periods_BK;


% compute moments, Om, as the autocovariances of the data for lags
% 0,1,...,K
K=4;
% Take the initial data, estimate a VAR
max_lags = 16;
% [AIC,BIC,HQ] = aic_bic_hq(filt_data',max_lags);
% p =min([AIC,BIC,HQ]);
% A is the impact matrix, identified via Cholesky, B is the beta_ols, res are
% the residuals, sigma is the estimated VC matrix.
[A,B,res,sigma] = sr_var(filt_data', p);

% Rewrite the VAR(p) as VAR(1) (see Hamilton, p. 273, Mac)
c = B(1,:); % coefficients of the constant
PHI = B(2:end,:)';
F = [PHI; [eye(nobs*(p-1)), zeros(nobs*(p-1),nobs)]];
Q = [[sigma, zeros(nobs,nobs*(p-1))]; zeros(nobs*(p-1),nobs*p)];
% check sizes
np = nobs*p;
sizeF = size(F) == [np,np];
sizeQ = size(Q) == [np,np];

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
Om = vec(Gamj);

res = (Om_data -Om).*diag(W1);
