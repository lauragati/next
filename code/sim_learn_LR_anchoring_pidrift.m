% simulate data from learning model LR learning of just constant
% with cemp's anchoring mechanism for the gain of PIBAR ONLY
% 18 sept 2019

function [xsim, ysim, shock, diff,pibar, k] = sim_learn_LR_anchoring_pidrift(gx,hx,eta,T,ndrop,e, Aa, Ab, As, param, setp,H, anal)
gbar = param.gbar;
ny = size(gx,1);
nx = size(hx,1);

ysim = zeros(ny,T);
xsim = zeros(nx,T);

%Learning PLM matrices, just a constant. Using RE as default starting point.
pibar = 0;
b = gx*hx;
brow_pibar = b(1,:);
bp = brow_pibar;

diff = zeros(T,1);
diff(1) = nan;
k = zeros(1,T);
k(:,1) = gbar^(-1);
%Simulate, with learning
for t = 1:T-1
    
    if t == 1
        ysim(:,t) = gx*xsim(:,t);
        xesim = hx*xsim(:,t);
    else
        %Form Expectations using last period's estimates
        if anal ==1
            [fa, fb] = fafb_anal_constant(param, setp, [pibar;0;0],b, xsim(:,t));
        else
            [fa, fb] = fafb_trunc_constant(param, setp, [pibar;0;0],b, xsim(:,t), H);
        end
        
        %Solve for current states
        ysim(:,t) = Aa*fa + Ab*fb + As*xsim(:,t);
        xesim = hx*xsim(:,t);
        
        %Update coefficients
        kt = fk_pidrift(pibar, b, xsim(:,t-1), k(:,t-1), param, setp, Aa, Ab, As);
        k(:,t) = kt;
        pibar = pibar + k(:,t).^(-1).*(ysim(1,t)-(pibar + bp*xsim(:,t-1)) );
        
        % check convergence
        diff(t) = max(max(abs(pibar - at_1)));
        
    end
    
    %Simulate transition with shock
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
