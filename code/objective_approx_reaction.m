function loss = objective_approx_reaction(coeffs,param,eN,capT,burnin,PLM,gain,gx,hx,eta,Aa,Ab,As, H, N)
bet  = param.bet;
alph = param.alph;
kapp = param.kapp;
lamx = param.lamx;

% Replace the estimated parameters with the new guess
param.('psi_pi') = coeffs(1);
param.('psi_x') = coeffs(2);

[ny,nx] = size(gx);

T = capT-H;
res = zeros(1,T);
% the innovations expected by the CB at each t
eN_CB = eN;
for t=2:T-1
    eN_CB(:,t+1:end,:) = zeros(nx,capT-t,N);
    
    % Simulate the model N times, at each t inputting the future innovations
    % the CB expects
    xsim = zeros(ny,capT,N);
    ysim = zeros(ny,capT,N);
    pibar = zeros(capT,N);
    k     = zeros(capT,N);
    g_pi  = zeros(capT,N);
    g_pibar = zeros(capT,N);
    
    
    for n=1:N
        e = squeeze(eN_CB(:,:,n));
        [xsim(:,:,n), ysim(:,:,n),phi_seq, k(:,n), g_pi(:,n), g_pibar(:,n)] = sim_learnLH_approx_reaction(gx,hx,eta,capT,burnin,e, Aa, Ab, As, param, PLM, gain);
        pibar(:,n) = squeeze(phi_seq(1,1,:));
    end
    
    b = gx*hx;
    
    % Take cross-sectional averages
    xsim_mean = mean(xsim,3);
    ysim_mean = mean(ysim,3);
    pibar_mean = mean(pibar,2);
    k_mean    = mean(k,2);
    g_pi_mean = mean(g_pi,2);
    g_pibar_mean = mean(g_pibar,2);
    
    % compute residual at t from target criterion
    pi = ysim_mean(1,:);
    x  = ysim_mean(2,:);
    pibar = pibar_mean;
    k = k_mean;
    g_pi  = g_pi_mean;
    g_pibar  = g_pibar_mean;
    % g_pi  = ones(size(g_pi)); a trick to see if things would be well-behaved
    % if g_pi wasn't blowing up
    s  = xsim_mean;
    
    % evaluate target criterion residual (comes from res_anchTC.m)
    sumprod = x(t+1); % when i=1
    for i=2:H-1
        prod = 1;
        for j=1:i-1
            prod = prod*(1-k(t+1+j)^(-1)-(pi(t+1+j)-pibar(t+j)-b(1,:)*s(:,t+j))*g_pibar(t+j)); % this assumed perfect foresight - not if I input the expected shocks in simulation
            %             prod = prod*(1-k(t+1+j)^(-1)-(pi(t+1+j)-pibar(t+j) -b(1,:)*hx^(j-1)*s(:,t))*g_pi(t+j)); % RE CB expectations for states
        end
        sumprod = sumprod + x(t+i)*prod;
    end
    
    res(t) = - pi(t) - lamx/kapp*x(t) ...
        +lamx/kapp*(1-alph)*bet/(1-alph*bet)*(k(t)^(-1) + (pi(t) - pibar(t-1) - b(1,:)*s(:,t-1))*g_pi(t))...
        *sumprod;
    
end

loss= max(abs(res));
