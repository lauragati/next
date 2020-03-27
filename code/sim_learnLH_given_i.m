% simulate data from long-horizon learning model given exog i-sequence
% based on sim_learnLH.m, but simplified:
% 1. took out any reference to endogenous states
% 2. took out any reference to forecasts and forecast errors
% 3. took out reference to whether expectations are anchored when shock hits
% 4. took out references to Jan 23-24 debate about how to write PLM-based forecasts
% 25 March 2020
function [pi,x,k,pibar,b,s,g_pi] = sim_learnLH_given_i(gx,hx,eta,capT,ndrop,e, Aa, Ab, As, param, PLM, gain, i_seq, dt, x0)

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

ysim = zeros(ny,capT);
xsim = zeros(nx,capT);

% Begin forcing ysim(3,:) to be the exogenous i-sequence
ysim(3,:) = i_seq;

if PLM==1 || PLM == 2 || PLM==11 || PLM == 21
    % 1 = Learning PLM matrices, just a constant.
    % 2 = learning slope and constant (HERE for entire vector of observables)
    % 11 = learning only a constant, just for inflation (scalar learning)
    % 21 = learning slope and constant, jsut for inflation (not implemented)
    a = zeros(ny,1);
    b = gx*hx; 
elseif PLM == -1
    % "only-mean" PLM
    % - ain't even using gx*hx for the forecast of inflation, only pibar
    a = zeros(ny,1);
    b = zeros(size(gx*hx));
else
    warning('I don''t know what you requested.')
end

learn_selector = [1,1,1]';
if PLM == 11 || PLM == 21
    learn_selector = [1,0,0]';
end
el = learn_selector;

% If endog gain, choose criterion
if gain == 21
    crit = 1; % CEMP's criterion
elseif gain == 22
    crit = 2; % CUSUM criterion
elseif gain == 23
    crit = 3; % smooth criterion (this is implemented ONLY for the case that only pi is learned)
end

phi = [a,b];
phi_seq = nan(ny,nx+1,capT);
phi_seq(:,:,1) = phi;

[sigy,sigx] = mom(gx,hx,eta*eta');
R = eye(nx+1); R(2:end,2:end) = sigx;
R_seq = repmat(R,[1,1,capT]);
R_seq(:,:,1) = R;

diff = zeros(capT,1);
diff(1) = nan;
k = zeros(1,capT);
k(:,1) = gbar^(-1);
g_pi = zeros(1,capT);


FA = nan(ny,capT);
FB = nan(ny,capT);

%%% initialize CUSUM variables: FEV om and criterion theta
% om = sigy; %eye(ny);
om = eta*eta';
% om = om(1,1);
thet = 0; % CEMP don't really help with this, but I think zero is ok.
% thet = thettilde; % actually it's quite sensitive to where you initialize it.
%%%

%Simulate, with learning
for t = 1:capT-1
    if t == 1
        ysim(:,t) = gx*xsim(:,t);
        xesim = hx*xsim(:,t);
    else
        % Select the variable that's being learned:
        a = el.*a;
        
        %Form Expectations using last period's estimates
        [fa, fb] = fafb_anal_constant_free(param,a, b, xsim(:,t),hx);
        FA(:,t) = fa; % save current LH expectations for output
        FB(:,t) = fb;
        
        %Solve for current states
        % this step needs to change!
        disp(['t=',num2str(t)])
%         dbstop in pi_x_given_i at 11 if t==33
        ysim(1:2,t) = pi_x_given_i(param,hx,fa,fb,xsim(:,t),i_seq(t));
%         ysim(:,t) = Aa*fa + Ab*fb + As*xsim(:,t);
        xesim = hx*xsim(:,t);
        
        %Update coefficients
        % Here the code differentiates between decreasing or constant gain
        if gain ==1 % decreasing gain
            k(:,t) = k(:,t-1)+1;
        elseif gain==21 || gain == 22 || gain == 23% endogenous gain
            if crit == 1 % CEMP's criterion
                fk = fk_CEMP(param,hx,Aa,Ab,As,a,b,eta,k(:,t-1));
            elseif crit==2 % CUSUM criterion
                fe = ysim(:,t)-(phi*[1;xsim(:,t-1)]); % short-run FE
%                 fe = fe(1,1);
                [fk, om, thet] = fk_CUSUM_vector(param,k(:,t-1),om, thet,fe);
%                 [fk, om, thet] = fk_cusum(param,k(:,t-1),om, thet,fe);
            elseif crit == 3 % smooth criterion
                fe = ysim(1,t)-(a(1) + b(1,:)*xsim(:,t-1)); 
                [fk, g_pi(t)] = fk_smooth_pi_only(param,fe,k(:,t-1));
            end
            k(:,t) = fk;
        elseif gain==3 % constant gain
            k(:,t) = gbar^(-1);
        end
        
        
        % Do the updating
        if PLM == 1 || PLM == -1 || PLM == 11 % when learning constant only, or "mean-only" PLM, or "constant-only, pi-only"
            a = el.*( a + k(:,t).^(-1).*( ysim(:,t)-(a + b*xsim(:,t-1)) )  );
            phi = [a,b];
        elseif PLM == 2 % constant and slope learning
            R = R + k(:,t).^(-1)*([1;xsim(:,t-1)]*[1;xsim(:,t-1)]' - R); % now I don't know if the gain should be the same here, Ryan uses the same gain.
            phi = (phi' + k(:,t).^(-1).*  (R\[1;xsim(:,t-1)] *(ysim(:,t)-phi*[1;xsim(:,t-1)])'))';
            
            % split phi into a and b
            a = phi(:,1);
            b = phi(:,2:end);
            
        end
        
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

% Split outputs
pi = ysim(1,:);
x  = xsim(1,:);
s = xsim;
pibar = squeeze(phi_seq(1,1,:));
