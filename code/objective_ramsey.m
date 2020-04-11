% the Ramsey problem's loss function
function loss = objective_ramsey(seq,param,gx,hx,eta,ndrop,eN)
bet = param.bet;
lamx = param.lamx;
lami = param.lami;

[ny,T,N] = size(eN);
ysim = zeros(ny,T,N);
for n=1:N
    e = squeeze(eN(:,:,n));
    [~, ysim(:,:,n)] = sim_learnLH_clean_smooth_given_seq(param,gx,hx,eta,seq,T,ndrop,e);
end

ysim2 = ysim.^2;
bet_t = zeros(T,1);
for t=1:T
    bet_t(t) = bet^t;
end
period_loss = bet_t.*( squeeze(ysim2(1,:,:) + lamx*ysim2(2,:,:) + lami*ysim2(3,:,:))  );
sum_loss = sum(period_loss,1);
expected_loss = mean(sum_loss);

loss = expected_loss;