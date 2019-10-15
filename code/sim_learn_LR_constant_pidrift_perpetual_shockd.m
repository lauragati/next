% same as sim_learn_LR_constant_pidrift_perpetual.m, except there's an innovation of
% x0 = vector of impulses with the innovation (delta)
% dt = innovation is imposed at time dt
% 14 oct 2019

function [xsim, ysim, shock, diff,pibar] = sim_learn_LR_constant_pidrift_perpetual_shockd(gx,hx,eta,T,ndrop,e, Aa, Ab, As, param, setp,H, anal,dt, x0)
gbar = param.gbar;

ny = size(gx,1);
nx = size(hx,1);

ysim = zeros(ny,T);
xsim = zeros(nx,T);

%Learning PLM matrices, just a constant. Using RE as default starting point.
pibar = 0;
b = gx*hx;
b1 = b(1,:);

diff = zeros(T,1);
diff(1) = nan;
k=gbar^(-1);
%Simulate, with learning
for t = 1:T-1
    
    if t == 1
        ysim(:,t) = gx*xsim(:,t);
        xesim = hx*xsim(:,t);
    else
        %Form Expectations using last period's estimates
        if anal ==1
            [fa, fb] = fafb_anal_constant(param, setp, [pibar;0;0], b, xsim(:,t));
        else
            [fa, fb] = fafb_trunc_constant(param, setp, [pibar;0;0], b, xsim(:,t), H);
        end
        
        %Solve for current states
        ysim(:,t) = Aa*fa + Ab*fb + As*xsim(:,t);
        xesim = hx*xsim(:,t);
        
        %Update coefficients
        pibar = pibar + k^(-1)* (ysim(1,t)-(pibar + b1*xsim(:,t-1)) );
        
        % check convergence
        diff(t) = max(max(abs(pibar - at_1)));
        
    end
    
    %Simulate transition with shock
    %%% here is the addition of the impulse
    if t+1==dt
        e(:,t+1) = e(:,t+1)+x0';
    end
    %%%
    xsim(:,t+1) = xesim + eta*e(:,t+1);
    
    % generate an old constant, to check convergence
    at_1 = pibar;
end

%Last period observables.
ysim(:,t+1) = gx*xsim(:,t+1);

%Drop ndrop periods from simulation
xsim = xsim(:,ndrop+1:end);
ysim = ysim(:,ndrop+1:end);
shock = e(:,ndrop+1:end); % innovations
