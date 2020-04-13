% the Ramsey problem's loss function
function loss = objective_ramsey_materials25(coeffs,param,gx,hx,eta,ndrop,eN)
bet = param.bet;
lamx = param.lamx;
lami = param.lami;

param.psi_pi = coeffs(1);
param.psi_x  = coeffs(2);
param.psi_k     = coeffs(3);
param.psi_pibar = coeffs(4);
param.psi_xbar  = coeffs(5);

% these are not really used, just too lazy to take 'em out
[Aa, Ab, As] = matrices_A_13_true_baseline(param, hx);

PLM = 1; % constant-only vector learning
[ny,T,N] = size(eN);
ysim = zeros(ny,T,N);
for n=1:N
    e = squeeze(eN(:,:,n));
    [~, ysim(:,:,n)] = sim_learnLH_clean_g(param,gx,hx,eta, Aa, Ab, As,PLM, T,ndrop,e);

end

ysim2 = ysim.^2;
bet_t = zeros(T,1);
for t=1:T
    bet_t(t) = bet^t;
end
period_loss = reshape(squeeze(ysim2(1,:,:) + lamx*ysim2(2,:,:) + lami*ysim2(3,:,:)),T,N) ;
discounted_period_loss = bet_t.* period_loss;
sum_loss = sum(discounted_period_loss,1);
expected_loss = mean(sum_loss);

loss = expected_loss;