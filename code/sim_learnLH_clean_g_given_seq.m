% sim_learnLH_clean_g
% vector learning for the smooth anchoring function g.
% also reverting to using inverse gains, k1, everywhere.
% based on sim_learnLH_clean.m
% 12 April 2020
function [xsim, ysim, k, phi_seq, FA, FB, diff] = sim_learnLH_clean_g_given_seq(param,gx,hx,eta, seq,PLM, T,ndrop,e, dt, x0)

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

learn_selector = [1,1,1]';
if PLM == 11 || PLM == 21
    learn_selector = [1,0,0]';
end
el = learn_selector;

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

%Simulate, with learning
for t = 1:T-1
    if t == 1
        ysim(:,t) = gx*xsim(:,t);
                ysim(end-size(seq,1)+1:end,t) = seq(:,t); % <-----
        xesim = hx*xsim(:,t);
    else
        % Select the variable that's being learned:
        a = el.*a;
        
        %Form Expectations using last period's estimates
        [fa, fb] = fafb_anal_constant_free(param,a, b, xsim(:,t),hx);
        FA(:,t) = fa;
        FB(:,t) = fb;
        
        %Solve for current states
        ysim(:,t) = A9A10(param,hx,fa,fb,xsim(:,t),seq(:,t)); % <-----
        xesim = hx*xsim(:,t);
        
        %Update coefficients
        % Only smooth criterion here
                % scalar:
                fe = ysim(1,t)-(a(1) + b(1,:)*xsim(:,t-1));
                fk = fk_smooth_pi_only(param,fe,k(:,t-1));
%                 % vector: 
%                 fe = ysim(:,t)-(phi*[1;xsim(:,t-1)]);
%                 fk = fk_smooth(param,fe,k(:,t-1));
            k(:,t) = fk;

        
        % Create forecasts and FE
        morning_fcst(:,t) = phi*[1;xsim(:,t)];
        % Yesterday evening's forecast error
        FEt_1(:,t) = ysim(:,t)-(phi*[1;xsim(:,t-1)]);
        
        % Do the updating
        if PLM == 1 || PLM == 11 % when learning constant only or "constant-only, pi-only"
            a = el.*( a + k(:,t).^(-1).*( ysim(:,t)-(a + b*xsim(:,t-1)) )  );
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
ysim(end-size(seq,1)+1:end,t+1) = seq(:,t+1); % <-----


%Drop ndrop periods from simulation
xsim = xsim(:,ndrop+1:end);
ysim = ysim(:,ndrop+1:end);
shock = e(:,ndrop+1:end); % innovations
k = k(:,ndrop+1:end);
