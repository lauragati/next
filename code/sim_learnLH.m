% simulate data from long-horizon learning model, no anchoring, but vector
% learning
% constant = 1 if learn constant only, -1 if "only-mean" learning, 2 if both constant and slope.
% gain: 1 = decreasing gain, 3 = constant gain
% dt = timing of shock. If 0 or not specified, then no shock.
% x0 = the shock.
% 18 Nov 2019

function [xsim, ysim, evening_fcst, morning_fcst, FA, FB, FEt_1, shock, diff,a_seq, k] = sim_learnLH(gx,hx,eta,T,ndrop,e, Aa, Ab, As, param, PLM, gain, dt, x0)

max_no_inputs = nargin('sim_learn');
if nargin < max_no_inputs %no shock specified
    dt = 0;
    x0 = 0;
end

gbar = param.gbar;
ny = size(gx,1);
nx = size(hx,1);

ysim = zeros(ny,T);
xsim = zeros(nx,T);

if PLM==1 
    %Learning PLM matrices, just a constant. Using RE as default starting point.
    a = zeros(ny,1);
    b = gx*hx;
elseif PLM == -1
    % "only-mean" PLM
    % - ain't even using gx*hx for the forecast of inflation, only pibar
    a = zeros(ny,1);
    b = zeros(size(gx*hx));
elseif PLM == 2
    % learning slope and constant (HERE for entire vector of observables)
    a = zeros(ny,1);
    b = gx*hx;
else
    warning('I don''t know what you requested.')
end

phi_seq = nan(T,nx+1);
a_seq = zeros(T,1);
diff = zeros(T,1);
diff(1) = nan;
k = zeros(1,T);
k(:,1) = gbar^(-1);
%%% initialize CUSUM variables: FEV om and criterion theta
om = 0;
thet = 0; % CEMP don't really help with this, but I think zero is ok.
%%%
evening_fcst = nan(T,1);
morning_fcst = nan(T,1);
FA = nan(ny,T);
FB = nan(ny,T);
FEt_1 = nan(T,1); % yesterday evening's forecast error, made at t-1 but realized at t and used to update pibar at t
%Simulate, with learning
for t = 1:T-1
    if t == 1
        ysim(:,t) = gx*xsim(:,t);
        xesim = hx*xsim(:,t);
    else
        %Form Expectations using last period's estimates
        [fa, fb] = fafb_anal_constant_free(param,[a;0;0], b, xsim(:,t),hx); % new hx version
        FA(:,t) = fa; % save current LH expectations for output
        FB(:,t) = fb;
        
        %Solve for current states
        ysim(:,t) = Aa*fa + Ab*fb + As*xsim(:,t);
        xesim = hx*xsim(:,t);
        
        %Update coefficients
        % Here the code differentiates between decreasing, endogenous or constant gain
        if gain ==1 % decreasing gain
            k(:,t) = k(:,t-1)+1;
        elseif gain==2 % endogenous gain
            % now choose the criterion
            if criterion == 1 % CEMP's
                kt = fk_pidrift_free([a;0;0], b, xsim(:,t-1), k(:,t-1), param, Aa, Ab, As, hx);  % new hx version
            elseif criterion == 2 % CUSUM
                % Cusum doesn't depend on P or n, so we need no difference
                % between free or not.
                fe = ysim(1,t)-(a + b1*xsim(:,t-1)); % short-run FE
                [kt, om, thet] = fk_cusum(param,k(:,t-1),omt_1, thett_1,fe);
            end
            k(:,t) = kt;
        elseif gain==3 % constant gain
            k(:,t) = gbar^(-1);
        end
        
        % Create forecasts, FE and do the updating
        if PLM == 1 || PLM == -1 % when learning constant only, or "mean-only" PLM
            morning_fcst(t) = a + b1*xsim(:,t); % this morning's one-step ahead forecast of tomorrow's state E(pi_{t+1} | I_{t}^m)
            FEt_1(t) = ysim(1,t)-(a + b1*xsim(:,t-1)); % yesterday evening's forecast error, realized today
            a = a + k(:,t).^(-1).*(ysim(1,t)-(a + b1*xsim(:,t-1)) );
            evening_fcst(t) = a + b1*xsim(:,t); % today's evening's one-step ahead forecast of tomorrow's state E(pi_{t+1} | I_{t}^e)
            
        elseif PLM == 2 % constant and slope learning
            phi = [a b1]; % gather constant and slope
            morning_fcst(t) = phi*[1 xsim(:,t)']'; 
            FEt_1(t) = ysim(1,t)-(phi*[1 xsim(:,t-1)']');
            phi = phi + k(:,t).^(-1).*(ysim(1,t)-(phi*[1 xsim(:,t-1)']')) ;
            evening_fcst(t) = phi*[1 xsim(:,t)']'; 
            % split phi into a and b once you've updated phi
            a = phi(1,1);
            b1 = phi(1,2:end);
            phi_seq(t,:) = phi; % store phis in case you wanna look at them
        end
        
        % check convergence
        diff(t) = max(max(abs(a - at_1)));
        
    end
    
    %Simulate transition with shock
    %%% here is the addition of the impulse
    if t+1==dt
        e(:,t+1) = e(:,t+1)+x0';
    end
    %%%
    xsim(:,t+1) = xesim + eta*e(:,t+1);
    
    % generate an old constant, to check convergence
    at_1 = a;
    a_seq(t)= a;
    %%% update CUSUM parameters
    omt_1 = om;
    thett_1 = thet;
    
end

%Last period observables.
ysim(:,t+1) = gx*xsim(:,t+1);

%Drop ndrop periods from simulation
xsim = xsim(:,ndrop+1:end);
ysim = ysim(:,ndrop+1:end);
shock = e(:,ndrop+1:end); % innovations
