% simulate data from learning model LR learning of both slope and constant
% Based on the structure of sim_learn_EE_check, except expectations are
% different.
% 15 sept 2019

function [xsim, ysim, shock] = sim_learn_LR(gx,hx,eta,T,ndrop,e, Aa, Ab, As, param, setp,H)

ny = size(gx,1);
nx = size(hx,1);

ysim = zeros(ny,T);
xsim = zeros(nx,T);

%Learning PLM matrices, include a constant. Using RE as default starting point. 
gl = [zeros(ny,1) gx*hx];
[~,sigx] = mom(gx,hx,eta*eta');
R = eye(nx+1); R(2:end,2:end) = sigx;

%Simulate, with learning
for t = 1:T-1
    
    if t == 1
        ysim(:,t) = gx*xsim(:,t);
        xesim = hx*xsim(:,t);
    else
        %Form Expectations using last period's estimates
        [fa_anal, fb_anal] = fafb_anal(param, setp, gl, xsim(:,t));
        [fa_trunc, fb_trunc] = fafb_trunc(param, setp, gl, xsim(:,t), H); 
        fa = fa_anal;
        fb = fb_anal;
        
        %Solve for current states
        ysim(:,t) = Aa*fa + Ab*fb + As*xsim(:,t);
        xesim = hx*xsim(:,t);
      
        %Update coefficients    
        R = R + t^(-1)*([1;xsim(:,t-1)]*[1;xsim(:,t-1)]' - R);
        gl = (gl' + t^(-1)*(R\([1;xsim(:,t-1)]*(ysim(:,t)-gl*[1;xsim(:,t-1)])')))';
        
    end
    
    %Simulate transition with shock
    xsim(:,t+1) = xesim + eta*e(:,t+1);
end

%Last period observables.
ysim(:,t+1) = gx*xsim(:,t+1);

%Drop ndrop periods from simulation
xsim = xsim(:,ndrop+1:end);
ysim = ysim(:,ndrop+1:end);
shock = e(:,ndrop+1:end); % innovations
