% a copy of sim_learnLH.m except it treats b as approximating gx, not gx*hx
% <--------- indicates the changes
% 25 Jan 2020

function [xsim, ysim, evening_fcst, morning_fcst, FA, FB, FEt_1, shock, diff,phi_seq, k] = sim_learnLH_materials14(gx,hx,eta,T,ndrop,e, Aa, Ab, As, param, PLM, gain, dt, x0)

this_code = mfilename;
max_no_inputs = nargin(this_code);
if nargin < max_no_inputs %no shock specified
    dt = 0;
    x0 = 0;
end

gbar = param.gbar;
ny = size(gx,1);
nx = size(hx,1);

ysim = zeros(ny,T);
xsim = zeros(nx,T);

if PLM==1 || PLM == 2
    % 1 = Learning PLM matrices, just a constant. Using RE as default starting point.
    % 2 = learning slope and constant (HERE for entire vector of observables)
    a = zeros(ny,1);
    %     b = gx*hx; % keep old formulation, just make sure you're consistent, 24 Jan 2020
    b = gx; % 23 Jan 2020 version <---------
elseif PLM == -1
    % "only-mean" PLM
    % - ain't even using gx*hx for the forecast of inflation, only pibar
    a = zeros(ny,1);
    b = zeros(size(gx*hx));
else
    warning('I don''t know what you requested.')
end

phi = [a,b];
phi_seq = nan(ny,nx+1,T);
phi_seq(:,:,1) = phi;

[~,sigx] = mom(gx,hx,eta*eta');
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
        
        % A mini-projection facility for slope-and-constant learning
        if PLM==2
            keep_R = projection_facility(R,0.995);
            if keep_R ~=1
                R = R_seq(:,:,t-1);
                phi = phi_seq(:,:,t-1);
            end
        end
        
        %Form Expectations using last period's estimates
        %         [fa, fb] = fafb_anal_constant_free(param,a, b, xsim(:,t),hx); % new hx version
        [fa, fb] = fafb_materials14(param,a, b, xsim(:,t),hx); % 23 Jan 2020 version <---------
        FA(:,t) = fa; % save current LH expectations for output
        FB(:,t) = fb;
        
        %Solve for current states
        ysim(:,t) = Aa*fa + Ab*fb + As*xsim(:,t);
        % Check if your other way of writing the ALM gives the same thing:
        [F,G] = FG_materials14(param,hx,Aa,Ab,As,a,b);
        ycheck = [F,G]*[1;xsim(:,t)];
        if max(abs(ycheck-ysim(:,t))) > 1e-12
            warning('ALM is not correct.')
        end
        xesim = hx*xsim(:,t);
        % If there are endogenous states...
        if endog_states==1
            %             % ...replace the last row of hx with the respective row of the estimated gx
            %             % xesim(4,1) = phi(lag_what,:)*[1;xsim(:,t)]; % 23 Jan 2020 version
            %             %             phx = pinv(hx);
            %             %             ghat = b*phx;
            %             %             ghat = inv(hx);
            %             %             xesim(4,1) = a(lag_what) + ghat(lag_what,:)*xsim(:,t); % 24 Jan 2020 version: need hx^-1
            %             xesim(4,1) = ysim(lag_what,t); % final version for now, 24 Jan 2020
            %             % HX = [[hx(1:3, 1:3), [0;0;0]] ; ghat(lag_what,:)];
            %             % max(abs(eig(HX)));
            %             % Projection facility following Ryan would be
            %             % keep_phi = projection_facility(HX,1);
%             HX = [[0;0;0;a(lag_what)], [hx(1:3, 1:4) ; b(lag_what,:)]]; % <---------
            HX = [[0;0;0;F(lag_what)], [hx(1:3, 1:4) ; G(lag_what,:)]]; % <--------- % <--------- NEW NEW thing
            xesim = HX*[1;xsim(:,t)]; % <---------
        end
        
        
        %Update coefficients
        % Here the code differentiates between decreasing or constant gain
        if gain ==1 % decreasing gain
            k(:,t) = k(:,t-1)+1;
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
            a = a + k(:,t).^(-1).*( ysim(:,t)-(a + b*hx*xsim(:,t-1)) ); % <---------
            phi = [a,b];
        elseif PLM == 2 % constant and slope learning
            R = R + k(:,t).^(-1)*([1;xsim(:,t-1)]*[1;xsim(:,t-1)]' - R); % now I don't know if the gain should be the same here, Ryan uses the same gain.
            phi = (phi' + k(:,t).^(-1).*  (R\[1;xsim(:,t-1)] *(ysim(:,t)-(a + b*hx*xsim(:,t-1)))'))'; % <---------
            
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
