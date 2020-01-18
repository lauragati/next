% simulate data from long-horizon learning model, no anchoring, but vector
% learning, but now where agents learn hx too!
% Use the same gain as for learning gx, and same R.
% constant = 1 if learn constant only, -1 if "only-mean" learning, 2 if both constant and slope.
% gain: 1 = decreasing gain, 3 = constant gain
% dt = timing of shock. If 0 or not specified, then no shock.
% x0 = the shock.
% 18 Jan 2020
function [xsim, ysim, evening_fcst, morning_fcst, FA, FB, FEt_1, shock, diff,phi_seq,diff_Gam,Gam_seq, k] = sim_learnLH_learnhx(gx,hx,eta,T,ndrop,e, Aa, Ab, As,Ba, Bb, param, PLM, gain, dt, x0)

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
% In any case initialize hx learning at RE: % NEW
c = zeros(nx,1);
d = hx;

phi = [a,b];
phi_seq = nan(ny,nx+1,T);
phi_seq(:,:,1) = phi;

Gam = [c,d]; % NEW
Gam_seq = nan(nx, nx+1,T); % NEW
Gam_seq(:,:,1) = Gam; % NEW

[~,sigx] = mom(gx,hx,eta*eta');
R = eye(nx+1); R(2:end,2:end) = sigx;

diff = zeros(T,1);
diff(1) = nan;
k = zeros(1,T);
k(:,1) = gbar^(-1);

evening_fcst = nan(ny,T);
morning_fcst = nan(ny,T);
FA = nan(ny,T);
FB = nan(ny,T);
GA = nan(nx,T);  % NEW
GB = nan(nx,T);  % NEW
FEt_1 = nan(ny,T); % yesterday evening's forecast error, made at t-1 but realized at t and used to update pibar at t
%Simulate, with learning
for t = 1:T-1
    if t == 1
        ysim(:,t) = gx*xsim(:,t);
        xesim = hx*xsim(:,t);
    else
        
        %Form Expectations using last period's estimates
        [fa, fb, ga, gb] = fafb_gagb(param, a, b, c, d, s, hx); % NEW
        FA(:,t) = fa; % save current LH expectations for output
        FB(:,t) = fb;
        GA(:,t) = ga;
        GB(:,t) = gb;
        
        %Solve for current states
        ysim(:,t) = Aa*fa + Ab*fb + As*xsim(:,t) + Ba*ga + Bb*gb;
        xesim = hx*xsim(:,t);
        
        %Update coefficients
        % Here the code differentiates between decreasing or constant gain
        if gain ==1 % decreasing gain
            k(:,t) = k(:,t-1)+1;
        elseif gain==3 % constant gain
            k(:,t) = gbar^(-1);
        end
        
        % Create forecasts, FE and do the updating
        if PLM == 1 || PLM == -1 % when learning constant only, or "mean-only" PLM
            morning_fcst(:,t) = phi*[1;xsim(:,t-1)];
            FEt_1(:,t) = ysim(:,t)-(phi*[1;xsim(:,t-1)]);
            a = a + k(:,t).^(-1).*( ysim(:,t)-(a + b*xsim(:,t-1)) );
            evening_fcst(:,t) = a + b*xsim(:,t-1);
            phi = [a,b];
            
            % update the hx-learning coefficients
            c = c + k(:,t).^(-1).*( xsim(:,t)-(c + d*xsim(:,t-1)) ); % NEW
            Gam = [c,d]; % NEW
        elseif PLM == 2 % constant and slope learning
            morning_fcst(:,t) = phi*[1;xsim(:,t-1)];
            FEt_1(:,t) = ysim(:,t)-(phi*[1;xsim(:,t-1)]);
            R = R + k(:,t).^(-1)*([1;xsim(:,t-1)]*[1;xsim(:,t-1)]' - R); % now I don't know if the gain should be the same here
            phi = (phi' + k(:,t).^(-1).*  (R\[1;xsim(:,t-1)] *(ysim(:,t)-phi*[1;xsim(:,t-1)])'))';
            evening_fcst(:,t) = phi*[1;xsim(:,t-1)];
            
            % split phi into a and b
            a = phi(:,1);
            b = phi(:,2:end);
            
            % update the hx-learning coefficients
            % using the same R since regressors are the same
            Gam = (Gam' + k(:,t).^(-1).*  (R\[1;xsim(:,t-1)] *(xsim(:,t)-Gam*[1;xsim(:,t-1)])'))';  % NEW
            
            % split into c and d
            c = Gam(:,1);
            d = Gam(:,2:end);
        end
        
        phi_seq(:,:,t) = phi; % store phis
        % check convergence
        diff(t) = max(max(abs(phi - squeeze(phi_seq(:,:,t-1)))));
        
         Gam_seq(:,:,t) = Gam; % NEW
         diff_Gam(t) = max(max(abs(Gam - squeeze(Gam_seq(:,:,t-1))))); % NEW
        
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