% sim_learnLH_clean_smooth.m
% a cleaned-up version of sim_learnLH.m, for the smooth (scalar) anchoring criterion
% only
% 10 April 2020
function [xsim, ysim, k, pibar, FA, FB, g_pi, g_pibar,fett_1eve, diff] = sim_learnLH_clean_smooth(param,gx,hx,eta, Aa, Ab, As, T,ndrop,e, dt, x0)

this_code = mfilename;
max_no_inputs = nargin(this_code);
if nargin < max_no_inputs %no shock specified
    dt = 0;
    x0 = 0;
else
    dt= dt+ndrop;
end

gbar = param.gbar;
ny = size(gx,1);
nx = size(hx,1);

ysim = zeros(ny,T);
xsim = zeros(nx,T);

a = zeros(ny,1);
b = gx*hx;
phi = [a,b];
phi_seq = nan(ny,nx+1,T);
phi_seq(:,:,1) = phi;

learn_selector = [1,0,0]';
el = learn_selector;

diff = zeros(T,1);
diff(1) = nan;
k = zeros(1,T);
k(:,1) = gbar^(-1);
pibar = zeros(1,T);
g_pi = zeros(1,T);
g_pibar = zeros(1,T);

FA = nan(ny,T);
FB = nan(ny,T);
fett_1eve = zeros(1,T);

%Simulate, with learning
for t = 1:T-1
    if t == 1
        ysim(:,t) = gx*xsim(:,t);
        xesim = hx*xsim(:,t);
    else
        % Select the variable that's being learned:
        a = el.*a;
        
        %Form Expectations using last period's estimates
        [fa, fb] = fafb_anal_constant_free(param,a, b, xsim(:,t),hx);
        FA(:,t) = fa;
        FB(:,t) = fb;
        
        %Solve for current states
        ysim(:,t) = Aa*fa + Ab*fb + As*xsim(:,t); 
%         ysim(:,t) = aleph_gimel(param,hx,fa,fb,xsim(:,t));
        xesim = hx*xsim(:,t);
        
        %Update coefficients
        fe = ysim(1,t)-(a(1) + b(1,:)*xsim(:,t-1));
        fett_1eve(t) = fe;
        [k(:,t), g_pi(:,t),g_pibar(:,t)] = fk_smooth_pi_only(param,fe,k(:,t-1));
        
        % Do the updating
        a = el.*( a + k(:,t).^(-1).*( ysim(:,t)-(a + b*xsim(:,t-1)) )  );
        phi = [a,b];
        
        phi_seq(:,:,t) = phi; % store phis
        pibar(t) = a(1,1); % store pibar-sequence
        % check convergence
        diff(t) = max(max(abs(phi - squeeze(phi_seq(:,:,t-1)))));       
    end
    
    %Simulate transition with shock
    %%% here is the addition of the impulse
    if t+1==dt
        e(:,t+1) = e(:,t+1)+x0';
    end
    %%%
    xsim(:,t+1) = xesim + eta*e(:,t+1);
    
end

%Last period observables.
ysim(:,t+1) = gx*xsim(:,t+1);

%Drop ndrop periods from simulation
xsim = xsim(:,ndrop+1:end);
ysim = ysim(:,ndrop+1:end);
shock = e(:,ndrop+1:end); % innovations
k = k(:,ndrop+1:end);
pibar = pibar(:,ndrop+1:end);
g_pi = g_pi(:,ndrop+1:end);
g_pibar = g_pibar(:,ndrop+1:end);
FA = FA(:,ndrop+1:end);
FB = FB(:,ndrop+1:end);

