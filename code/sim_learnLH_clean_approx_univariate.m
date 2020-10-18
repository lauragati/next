% a version of sim_learnLH_clean.m where the additional input is an
% approximation coefficients matrix (alph) to approximate the
% evolution of the gain using Ryan's ndim_simplex (ND piecewise linear interpolation)
% In addition, this version does that for the univariate anchoring function
% 21 June 2020
% Udpdate 22 July 2020: added measurement error v
function [xsim, ysim, k, phi_seq, FA, FB, FEt_1,diff, explode_count, negk_count, explode_t, negk_t, k_dt] = sim_learnLH_clean_approx_univariate(alph,x,param,gx,hx,eta, PLM, gain, T,ndrop,e,v,...
    knowTR,mpshock, dt, x0)

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

% If endog gain, choose criterion
if gain == 21
    crit = 1; % CEMP's criterion
elseif gain == 22
    crit = 2; % CUSUM criterion
elseif gain == 23
    crit = 3; % smooth criterion 
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
FA = zeros(ny,T);
FB = zeros(ny,T);
FEt_1 = zeros(ny,T); % yesterday evening's forecast error, made at t-1 but realized at t and used to update pibar at t

%%% initialize CUSUM variables: FEV om and criterion theta
% om = sigy; %eye(ny);
om = eta*eta';
% om = om(1,1);
thet = 0; % CEMP don't really help with this, but I think zero is ok.
% thet = thettilde; % actually it's quite sensitive to where you initialize it.
%%%
explode_count=0;
negk_count =0;
explode_t = zeros(T,1);
negk_t = zeros(T,1);


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
        ysim(:,t) = ALM(param,hx,fa,fb,xsim(:,t), knowTR,mpshock);
        xesim = hx*xsim(:,t);
        
        %Update coefficients
        % Here the code differentiates between decreasing or constant gain
        if gain ==1 % decreasing gain
            k(:,t) = k(:,t-1)+1;
        elseif gain==21 || gain == 22 || gain == 23% endogenous gain
            if crit == 1 % CEMP's criterion
                fk = fk_CEMP(param,hx,a,b,eta,k(:,t-1));
            elseif crit==2 % CUSUM criterion
                fe = ysim(:,t)-(phi*[1;xsim(:,t-1)]); % short-run FE
                %                 fe = fe(1,1);
                [fk, om, thet] = fk_CUSUM_vector(param,k(:,t-1),om, thet,fe);
                %                 [fk, om, thet] = fk_cusum(param,k(:,t-1),om, thet,fe);
            elseif crit == 3 % smooth criterion
                fe = ysim(1,t)-(a(1) + b(1,:)*xsim(:,t-1));
%                 fex = ysim(2,t)-(a(2) + b(2,:)*xsim(:,t-1)); % output gap fe

                % let's try to catch explosive paths
                fe_thresh = 5;
                if abs(fe) > fe_thresh
                    tfek1 = [t;fe; inv(k(:,t-1))];
%                     if explode_count == 0
%                         disp(['Setting fe to upper threshold. t, fe, k1_{t-1} : ', num2str(tfek1'), '; and obs: ', num2str(ysim(:,t-1)')])
%                     end
                    fe = sign(fe)*fe_thresh;
                    explode_count=explode_count+1;
                    explode_t(t) =t;
                end
                
                fk = fk_smooth_approx_univariate(alph,x,fe);
                

                
%                 fk_thresh = 5.8824; % corresponds to a gain of 0.17
%                 if fk < fk_thresh
%                     tfek1 = [t;fe; inv(fk)];
%                     if already_exploded == 0
%                         disp(['Setting fk to a lower threshold. t, fe, k1_t : ', num2str(tfek1'), '; and obs: ', num2str(ysim(:,t-1)')])
%                     end
%                     fk = fk_thresh;
%                     already_exploded=1;
%                 end

                
            end
            k(:,t) = fk;
            if k(:,t) < 0
                disp(['Negative gain k^-1 = ', num2str(inv(k(:,t))), ', fe = ', num2str(fe)])
                disp(['t = ', num2str(t)])
                disp(['k^-1_{t-1} = ', num2str(inv(k(:,t-1)))])
                disp('setting k^-1 to 0')
                k(:,t) = 10^16;
                negk_count = negk_count+1;
                negk_t(t) =t;

            end
        elseif gain==3 % constant gain
            k(:,t) = gbar^(-1);
%             k(:,t) = 0.0478^(-1); % to do IRFs with cgain learning with
%             gain at simulated median
        end
        
        % Create forecasts and FE
        morning_fcst(:,t) = phi*[1;xsim(:,t)];
        % Yesterday evening's forecast error
        FEt_1(:,t) = ysim(:,t)-(phi*[1;xsim(:,t-1)]);
        
        % Do the updating
        if PLM == 1 || PLM == 11 % when learning constant only or "constant-only, pi-only"
            a = el.*( a + k(:,t).^(-1).*( ysim(:,t)-(a + b*xsim(:,t-1)) )); % default
            % try fe_{t|t-2} here
%             if t>2
%             a = el.*( a + k(:,t).^(-1).*( ysim(:,t)-(squeeze(phi_seq(1,1,t-2)) + b*xsim(:,t-1)) )  ); 
%             else
%                 a = el.*( a + k(:,t).^(-1).*( ysim(:,t)-(a + b*xsim(:,t-1)) )  );
%             end
            phi = [a,b];
        elseif PLM == 2 % constant and slope learning
            R = R + k(:,t).^(-1)*([1;xsim(:,t-1)]*[1;xsim(:,t-1)]' - R); % now I don't know if the gain should be the same here, Ryan uses the same gain.
            phi = (phi' + k(:,t).^(-1).*  (R\[1;xsim(:,t-1)] *(ysim(:,t)-phi*[1;xsim(:,t-1)])'))';
            
            % split phi into a and b
            a = phi(:,1);
            b = phi(:,2:end);
            
        end
        evening_fcst(:,t) = phi*[1;xsim(:,t)]; % 23 Jan 2020 version
        
        phi_seq(:,:,t) = phi; % store phis 
        % check convergence
        diff(t) = max(max(abs(phi - squeeze(phi_seq(:,:,t-1)))));
        
    end
    
    %Simulate transition with shock
    %%% here is the addition of the impulse
    if t+1==dt % dt = dt+ndrop (see top of code)
        e(:,t+1) = e(:,t+1)+x0';
    elseif t==dt
        k_dt = k(:,t); % gather the gain when the shock hits 
    end
    %%%
    xsim(:,t+1) = xesim + eta*e(:,t+1);
    
end

%Last period observables.
ysim(:,t+1) = gx*xsim(:,t+1);

% Add measurement error to observables
ysim = ysim + v(1:3,:);
if size(v,1)>ny
phi_seq(1,1,:) = squeeze(phi_seq(1,1,:)) + v(4,:)';
end

%Drop ndrop periods from simulation
xsim = xsim(:,ndrop+1:end);
ysim = ysim(:,ndrop+1:end);
shock = e(:,ndrop+1:end); % innovations
k = k(:,ndrop+1:end);
FEt_1 = FEt_1(:,ndrop+1:end);
phi_seq = phi_seq(:,:,ndrop+1:end);
