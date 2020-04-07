% Optimizing i_seq. Tricky with the time dimension b/c
% if we choose T i-values optimally, for the last one we need the T...T+H
% future values, where H itself is a truncation (ideally it should be infty)
% Inputs
% i_seq = the sequence of i(t), t=1,...,T, T+1, ... T+H we're optimizing over
% (for now the t=T+1...T+H are only semi-optimal because we compute the target criterion only for t=1,...,T)
% 24 March 2020

function loss = objective_target_criterion(iseq,param,e,T,burnin,PLM,gain, gx,hx,SIG,Aa,Ab,As)
lamx = param.lamx;
kapp = param.kapp;
bet  = param.bet;
alph = param.alph;

capT = length(iseq);
H = capT - T;
% pi,x,k, pibar,s, g_pi are (T+H x 1)
% [pi,x,k,pibar,b,s,g_pi] = fun_sim_anchoring_given_i(param,capT,N,burnin,e,PLM,gain,iseq);
[s, y, ~, ~, ~, ~, ~, ~, ~, phi_seq, k, ~, g_pi] = sim_learnLH_given_i(gx,hx,SIG,capT+burnin,burnin,e, Aa, Ab, As, param, PLM, gain, iseq);
pi = y(1,:);
x  = y(2,:);
b = gx*hx;
pibar = squeeze(phi_seq(:,1,:));

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

% loss = resid'*resid; % squared residuals
loss = max(abs(resid)); % max norm of the residual
