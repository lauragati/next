% simulate data from long-horizon learning model
% general code, with lots of subcases
% constant = 1 if learn constant only - right now just give warning if you
% input something else
% gain: 1 = decreasing gain, 2 = endogenous gain, 3 = constant gain
% criterion: 1= CEMP's, 2 = CUSUM
% dt = timing of shock. If 0 or not specified, then no shock.
% x0 = the shock.
% free=1; % use versions of the code that are n-free (use hx instead of P)
% anchored: % 1 denotes anchored, 0 denotes unanchored when shock hits
% 19 Oct 2019

function [xsim, ysim, evening_fcst, morning_fcst, FA, FB, shock, diff,pibar_seq, k, anchored_when_shock] = sim_learn(gx,hx,eta,T,ndrop,e, Aa, Ab, As, param, setp,H, anal, constant, gain, criterion, free, dt, x0)
if nargin < 18 %no shock specified
    dt = 0;
    x0 = 0;
end
if nargin < 17
    free=1; % make using hx the default case instead of P.
end
anchored_when_shock = nan; % in case you call it w/ exogenous gain
gbar = param.gbar;
ny = size(gx,1);
nx = size(hx,1);

ysim = zeros(ny,T);
xsim = zeros(nx,T);

if constant==1
    %Learning PLM matrices, just a constant. Using RE as default starting point.
    pibar = 0;
    b = gx*hx;
    b1 = b(1,:);
else
    warning('Requested slope learning, but only learning constant is implemented.')
end

pibar_seq = zeros(T,1);
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
%Simulate, with learning
for t = 1:T-1
    
    if t == 1
        ysim(:,t) = gx*xsim(:,t);
        xesim = hx*xsim(:,t);
    else
        %Form Expectations using last period's estimates
        if anal ==1
            if free==0
                [fa, fb] = fafb_anal_constant(param, setp, [pibar;0;0], b, xsim(:,t)); %% old version with P
            elseif free==1
                [fa, fb] = fafb_anal_constant_free(param, setp, [pibar;0;0], b, xsim(:,t),hx); % new hx version
                FA(:,t) = fa; % save current LH expectations for output
                FB(:,t) = fb;
            end
        else
            [fa, fb] = fafb_trunc_constant(param, setp, [pibar;0;0], b, xsim(:,t), H);
        end
        
        %Solve for current states
        ysim(:,t) = Aa*fa + Ab*fb + As*xsim(:,t);
        xesim = hx*xsim(:,t);
        
        %Update coefficients
        % Here the code differentiates between decreasing, endogenous or
        % constant gain
        if gain ==1
            k(:,t) = k(:,t-1)+1;
        elseif gain==2
            % now choose the criterion
            if criterion == 1 % CEMP's
                if free==0
                    kt = fk_pidrift([pibar;0;0], b, xsim(:,t-1), k(:,t-1), param, setp, Aa, Ab, As); %% old version with P
                elseif free==1
                    kt = fk_pidrift_free([pibar;0;0], b, xsim(:,t-1), k(:,t-1), param, setp, Aa, Ab, As, hx);  % new hx version
                end
            elseif criterion == 2 % CUSUM
                % Cusum doesn't depend on P or n, so we need no difference
                % between free or not.
                f = ysim(1,t)-(pibar + b1*xsim(:,t-1)); % short-run FE
                [kt, om, thet] = fk_cusum(param,k(:,t-1),omt_1, thett_1,f);
            end
            k(:,t) = kt;
        elseif gain==3
            k(:,t) = gbar^(-1);
        end
        
        morning_fcst(t) = pibar + b1*xsim(:,t); % this morning's one-step ahead forecast of tomorrow's state E(pi_{t+1} | I_{t}^m)
        pibar = pibar + k(:,t).^(-1).*(ysim(1,t)-(pibar + b1*xsim(:,t-1)) );
        evening_fcst(t) = pibar + b1*xsim(:,t); % today's evening's one-step ahead forecast of tomorrow's state E(pi_{t+1} | I_{t}^e)
        % check convergence
        diff(t) = max(max(abs(pibar - at_1)));
        
    end
    
    %Simulate transition with shock
    %%% here is the addition of the impulse
    if t+1==dt
%         disp(['shock imposed at ', num2str(dt)])
        % check if anchored or not when shock hits
        if k(:,t) == gbar^(-1)
            anchored_when_shock = 0;
        elseif k(:,t) > gbar^(-1)
            anchored_when_shock = 1;
        end
        e(:,t+1) = e(:,t+1)+x0';
    end
    %%%
    xsim(:,t+1) = xesim + eta*e(:,t+1);
    
    % generate an old constant, to check convergence
    at_1 = pibar;
    pibar_seq(t)= pibar;
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
