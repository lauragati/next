function EL = loss_CB(param,y)
% expected loss of CB given a cross section of simulated histories
% y = (ny,T,N)
% 26 Jan 2020
lamx = param.lamx;
lami = param.lami;
bet = param.bet;

[~,T,N] = size(y);
pi_seq = squeeze(y(1,:,:));
x_seq = squeeze(y(2,:,:));
i_seq = squeeze(y(3,:,:));

X2 = x_seq.^2;
PI2 = pi_seq.^2;
I2  = i_seq.^2;

loss = zeros(1,N);
for n=1:N
    period_loss = zeros(1,T);
    for t=1:T
        period_loss(t) = bet^t * (lamx*X2(t,n) + PI2(t,n) + lami*I2(t,n));
    end
    loss(n) = 1/2*sum(period_loss);
end
EL = mean(loss);