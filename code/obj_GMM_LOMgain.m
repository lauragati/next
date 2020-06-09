function loss = obj_GMM_LOMgain(coeffs,param,e,T,ndrop,PLM,gain, Om_data, W1)

% Replace the estimated parameters with the new guess
pp.coefs = coeffs;

sig_r = param.sig_r;
sig_i = param.sig_i;
sig_u = param.sig_u;

% RE model
[fyn, fxn, fypn, fxpn] = model_NK(param);
[gx,hx]=gx_hx_alt(fyn,fxn,fypn,fxpn);
[ny, nx] = size(gx);
[Aa, Ab, As] = matrices_A_13_true_baseline(param, hx);
SIG = eye(nx).*[sig_r, sig_i, sig_u]';
eta = SIG; %just so you know

% Simulate data given parameters
[~, y] = sim_learnLH_clean(param,gx,hx,eta, PLM, gain, T,ndrop,e);

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

% Compute the model-implied moments
% compute moments, Om, as the autocovariances of the data for lags
% 0,1,...,K
K=4;
% Take the initial data, estimate a VAR
max_lags = 16;
[AIC,BIC,HQ] = aic_bic_hq(filt',max_lags);
p =min([AIC,BIC,HQ]);
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