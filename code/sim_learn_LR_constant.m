% simulate data from learning model LR learning of just constant
% 17 sept 2019

function [xsim, ysim, shock, diff,a] = sim_learn_LR_constant(gx,hx,eta,T,ndrop,e, Aa, Ab, As, param, setp,H, anal)

ny = size(gx,1);
nx = size(hx,1);

ysim = zeros(ny,T);
xsim = zeros(nx,T);

%Learning PLM matrices, just a constant. Using RE as default starting point. 
a = zeros(ny,1);
b = gx*hx;

diff = zeros(T,1);
diff(1) = nan;
%Simulate, with learning
for t = 1:T-1
    
    if t == 1
        ysim(:,t) = gx*xsim(:,t);
        xesim = hx*xsim(:,t);
    else
        %Form Expectations using last period's estimates
        if anal ==1
        [fa, fb] = fafb_anal_constant(param, setp, a,b, xsim(:,t));
        else
        [fa, fb] = fafb_trunc_constant(param, setp, a,b, xsim(:,t), H); 
        end
        
        %Solve for current states
        ysim(:,t) = Aa*fa + Ab*fb + As*xsim(:,t);
        xesim = hx*xsim(:,t);
      
         %Update coefficients    
        a = a + t^(-1)* (ysim(:,t)-(a + b*xsim(:,t-1)) );
        
        % check convergence
        diff(t) = max(max(abs(a - at_1)));
        
    end
    
    %Simulate transition with shock
    xsim(:,t+1) = xesim + eta*e(:,t+1);
    
    % generate an old constant, to check convergence
    at_1 = a;
end

%Last period observables.
ysim(:,t+1) = gx*xsim(:,t+1);

%Drop ndrop periods from simulation
xsim = xsim(:,ndrop+1:end);
ysim = ysim(:,ndrop+1:end);
shock = e(:,ndrop+1:end); % innovations