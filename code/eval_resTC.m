function resid = eval_resTC(param,ysim,k,phi_seq,s,g_pi,H)
lamx = param.lamx;
kapp = param.kapp;
bet  = param.bet;
alph = param.alph;

T = length(k)-H;
pi = ysim(1,:);
x  = ysim(2,:);
pibar = squeeze(phi_seq(:,1,:));
b = squeeze(phi_seq(:,2:end,1));

resid = zeros(T,1);
for t=2:T
sumprod = x(t+1);
for i=2:H
    prod = 1;
    for j=1:i-1
        prod = prod*(1-k(t+1+j)^(-1)-(pi(t+1+j)-pibar(t+j) -b(1,:)*s(:,t+j))*g_pi(t+j));
    end
    sumprod = sumprod + x(t+i)*prod;
end

resid(t) = - pi(t) - lamx/kapp*x(t) ...
    +lamx/kapp*(1-alph)*bet/(1-alph*bet)*(k(t)^(-1) + (pi(t) - pibar(t-1) - b(1,:)*s(:,t-1))*g_pi(t))...
    *sumprod;
end