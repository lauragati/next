% separate, hopefully correct, objective function to evaluate the residuals
% from the anchoring target criterion
% 9 April 2020
function loss = objective_seq_anchTC(seq,param,e,capT,burnin,PLM,gain,gx,hx,eta,Aa,Ab,As, H)
lamx = param.lamx;
kapp = param.kapp;
bet  = param.bet;
alph = param.alph;
sig  = param.sig;
nx = length(hx);

T = capT-H;
[ny,nx] = size(gx);
b = gx*hx;

stuff1 = [sig,1-bet,-sig*bet]; % fa(3) I hope is rational
% stuff1 = [sig,1-bet,0];
stuff2 = sig*[1,0,0]*(eye(nx)-bet*hx)^(-1);
stuff3 = [(1-alph)*bet,kapp*alph*bet,0];
stuff4 = [0,0,1]*(eye(nx)-alph*bet*hx)^(-1);

resTC = zeros(T,1);
resA9 = zeros(T,1);
resA10 = zeros(T,1);

% the innovations expected by the CB at each t
e_CB = e;
for t=2:T-1
    e_CB(:,t+1:end,:) = zeros(nx,capT-t);
    
    % Simulate the model 1 time, at each t inputting the future innovations
[xsim, ysim, ~, ~, fa, fb, ~, ~, ~,phi_seq, k,~,g_pi, g_pibar]= sim_learnLH_given_seq(gx,hx,eta,capT,burnin,e, Aa, Ab, As, param, PLM, gain,seq, H, 3);
    pibar = squeeze(phi_seq(1,1,:));
    
    % compute residual at t from target criterion
    pi = ysim(1,:);
    x  = ysim(2,:);
    s  = xsim;
    
    % 1.) evaluate target criterion residual (comes from res_anchTC.m)
    sumprod = x(t+1); % when i=1
    for i=2:H-1
        prod = 1;
        for j=1:i-1
            prod = prod*(1-k(t+1+j)^(-1)-(pi(t+1+j)-pibar(t+j)-b(1,:)*s(:,t+j))*g_pibar(t+j)); % this assumed perfect foresight - not if I input the expected shocks in simulation
            %             prod = prod*(1-k(t+1+j)^(-1)-(pi(t+1+j)-pibar(t+j) -b(1,:)*hx^(j-1)*s(:,t))*g_pi(t+j)); % RE CB expectations for states
        end
        sumprod = sumprod + x(t+i)*prod;
    end
    
    resTC(t) = - pi(t) - lamx/kapp*x(t) ...
        +lamx/kapp*(1-alph)*bet/(1-alph*bet)*(k(t)^(-1) + (pi(t) - pibar(t-1) - b(1,:)*s(:,t-1))*g_pi(t))...
        *sumprod;
    % 2.) Evaluate IS and PC residuals
    resA9(t) = -x(t) -sig*ysim(3,t) + stuff1*fb(:,t) + stuff2*s(:,t);
    resA10(t) =-pi(t) +kapp*x(t) + stuff3*fa(:,t) + stuff4*s(:,t);
    
end
resids = [resA9'; resA10'; resTC'];

loss= max(max(abs(resids)));