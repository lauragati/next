% simulate data from long-horizon learning model, no anchoring, but vector
% learning
% constant = 1 if learn constant only, -1 if "only-mean" learning, 2 if both constant and slope.
% gain: 1 = decreasing gain, 3 = constant gain
% dt = timing of shock. If 0 or not specified, then no shock.
% x0 = the shock.
% 18 Nov 2019

function [xsim, ysim, evening_fcst, morning_fcst, FA, FB, FEt_1, shock, diff,phi_seq, k] = sim_learnLH(gx,hx,eta,T,ndrop,e, Aa, Ab, As, param, PLM, gain, dt, x0)

this_code = mfilename;
max_no_inputs = nargin(this_code);
if nargin < max_no_inputs %no shock specified
    dt = 0;
    x0 = 0;
else
    dt= dt+ndrop;
end

gbar = param.gbar;
thettilde = param.thettilde;
ny = size(gx,1);
nx = size(hx,1);

ysim = zeros(ny,T);
xsim = zeros(nx,T);

if PLM==1 || PLM == 2
    % 1 = Learning PLM matrices, just a constant. Using RE as default starting point.
    % 2 = learning slope and constant (HERE for entire vector of observables)
    a = zeros(ny,1);
    b = gx*hx; % keep old formulation, just make sure you're consistent, 24 Jan 2020
    %     b = gx; % 23 Jan 2020 version
elseif PLM == -1
    % "only-mean" PLM
    % - ain't even using gx*hx for the forecast of inflation, only pibar
    a = zeros(ny,1);
    b = zeros(size(gx*hx));
else
    warning('I don''t know what you requested.')
end

% If endog gain, choose criterion
if gain == 21
    crit = 1; % CEMP's criterion
elseif gain == 22
    crit = 2; % CUSUM criterion
end

phi = [a,b];
phi_seq = nan(ny,nx+1,T);
phi_seq(:,:,1) = phi;

[sigy,sigx] = mom(gx,hx,eta*eta');
R = eye(nx+1); R(2:end,2:end) = sigx;
R_seq = repmat(R,[1,1,T]);
R_seq(:,:,1) = R;

diff = zeros(T,1);
diff(1) = nan;
k = zeros(1,T);
k(:,1) = gbar^(-1);

evening_fcst = nan(ny,T);
morning_fcst = nan(ny,T);
FA = nan(ny,T);
FB = nan(ny,T);
FEt_1 = nan(ny,T); % yesterday evening's forecast error, made at t-1 but realized at t and used to update pibar at t

%%% initialize CUSUM variables: FEV om and criterion theta
% om = sigy; %eye(ny);
om = eta*eta';
% om = om(1,1);
thet = 0; % CEMP don't really help with this, but I think zero is ok.
% thet = thettilde; % actually it's quite sensitive to where you initialize it.
%%%

% Do an initial check to see whether we have endogenous states that are
% lagged jumps
% If there are endogenous states...
if nx == 4
    endog_states=1;
    % ... check which ones
    lag_what = 0;
    for i=1:ny
        if max(abs((hx(4,:) - gx(i,:)))) < 1e-14
            lag_what = i;
        end
    end
    if lag_what==0
        warning('Couldn''t identify the lagged jump.')
    end
else
    endog_states=0;
end

%Simulate, with learning
for t = 1:T-1
    if t == 1
        ysim(:,t) = gx*xsim(:,t);
        xesim = hx*xsim(:,t);
    else
        
%         % A mini-projection facility for slope-and-constant learning
%         if PLM==2
%             keep_R = projection_facility(R,0.995);
%             if keep_R ~=1
%                 R = R_seq(:,:,t-1);
%                 phi = phi_seq(:,:,t-1);
%                 warning(['Projection facility invoked at t=', num2str(t) ])
%             end
%         end
        
        %Form Expectations using last period's estimates
        [fa, fb] = fafb_anal_constant_free(param,a, b, xsim(:,t),hx); % new hx version
        %         [fa, fb] = fafb_materials14(param,a, b, xsim(:,t),hx); % 23 Jan 2020 version
        FA(:,t) = fa; % save current LH expectations for output
        FB(:,t) = fb;
        
        %Solve for current states
        ysim(:,t) = Aa*fa + Ab*fb + As*xsim(:,t);
        xesim = hx*xsim(:,t);
        % If there are endogenous states...
        if endog_states==1
            % ...replace the last row of hx with the respective row of the estimated gx
            % xesim(4,1) = phi(lag_what,:)*[1;xsim(:,t)]; % 23 Jan 2020 version
            phx = pinv(hx);
            ghat = b*phx;
%             ghat = inv(hx);
            xesim(4,1) = a(lag_what) + ghat(lag_what,:)*xsim(:,t); % 24 Jan 2020 version: I've checked and this works (corresponds to replacing xsim(e,t+1) with ysim(e,t))
        end
        
        
        %Update coefficients
        % Here the code differentiates between decreasing or constant gain
        if gain ==1 % decreasing gain
            k(:,t) = k(:,t-1)+1;
        elseif gain==21 || gain == 22 % endogenous gain
            if crit == 1 % CEMP's criterion
                fk = fk_CEMP(param,hx,Aa,Ab,As,a,b,eta,k(:,t-1));
            elseif crit==2 % CUSUM criterion
                fe = ysim(:,t)-(phi*[1;xsim(:,t-1)]); % short-run FE
%                 fe = fe(1,1);
                [fk, om, thet] = fk_CUSUM_vector(param,k(:,t-1),om, thet,fe);
%                 [fk, om, thet] = fk_cusum(param,k(:,t-1),om, thet,fe);
            end
            k(:,t) = fk;
        elseif gain==3 % constant gain
            k(:,t) = gbar^(-1);
        end
        
        % Create forecasts and FE
        morning_fcst(:,t) = phi*[1;xsim(:,t)]; % 23 Jan 2020 version
        % morning_fcst(:,t) = phi*[1;xsim(:,t-1)]; % this should be xsim(:,t)
        % Yesterday evening's forecast error
        FEt_1(:,t) = ysim(:,t)-(phi*[1;xsim(:,t-1)]);
        
        % Do the updating
        if PLM == 1 || PLM == -1 % when learning constant only, or "mean-only" PLM
            a = a + k(:,t).^(-1).*( ysim(:,t)-(a + b*xsim(:,t-1)) );
            phi = [a,b];
        elseif PLM == 2 % constant and slope learning
            R = R + k(:,t).^(-1)*([1;xsim(:,t-1)]*[1;xsim(:,t-1)]' - R); % now I don't know if the gain should be the same here, Ryan uses the same gain.
            phi = (phi' + k(:,t).^(-1).*  (R\[1;xsim(:,t-1)] *(ysim(:,t)-phi*[1;xsim(:,t-1)])'))';
            
            % split phi into a and b
            a = phi(:,1);
            b = phi(:,2:end);
            
        end
        evening_fcst(:,t) = phi*[1;xsim(:,t)]; % 23 Jan 2020 version
        % evening_fcst(:,t) = phi*[1;xsim(:,t-1)]; % this should be xsim(:,t)
        
        phi_seq(:,:,t) = phi; % store phis
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
