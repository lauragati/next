% simulate data from learning model with Euler equation learning
% Rely strongly on Ryan's code sim_learn.m - here agents learn constant and
% slope
% 14 sept 2019

function [xsim, ysim, shock] = sim_learn_EE(gx,hx,fxp,fx,fyp,fy,eta,T,ndrop,e,param,PLM, gain, dt,x0)

max_no_inputs = nargin('sim_learn_EE');
if nargin < max_no_inputs %no shock specified
    dt = 0;
    x0 = 0;
end

gbar = param.gbar;
ny = size(gx,1);
nx = size(hx,1);

ysim = zeros(ny,T);
xsim = zeros(nx,T);

%Learning PLM matrices, include a constant. Using RE as default starting point.
gl = [zeros(ny,1) gx*hx];
[~,sigx] = mom(gx,hx,eta*eta');
R = eye(nx+1); R(2:end,2:end) = sigx;

a = zeros(ny,1);
b = gx*hx;

k = zeros(1,T);
k(:,1) = gbar^(-1);

detR=4*ones(T-1,1);
detR(1) = det(R);
gl_seq = gl;
R_seq  = R;
fake_eigs = eig((gl*gl').^(1/2));
fake_eigs_seq = fake_eigs;
max_mod = max(abs(fake_eigs));
max_mod_seq = max_mod;
%Simulate, with learning
for t = 1:T-1
    
    if t == 1
        ysim(:,t) = gx*xsim(:,t);
        xesim = hx*xsim(:,t);
    else
        %Form Expectations using last period's estimates
        yp_e = gl*[1; xsim(:,t)];
        yx = -[fy fxp]\(fyp*yp_e + fx*xsim(:,t));
        
        %Solve for current states
        ysim(:,t) = yx(1:ny);
        xesim     = yx(ny+1:end);
        
        %Update coefficients
        % Here the code differentiates between decreasing or constant gain
        if gain ==1 % decreasing gain
            k(:,t) = k(:,t-1)+1;
        elseif gain==3 % constant gain
            k(:,t) = gbar^(-1);
        end
        
        detR(t) = det(R);
        if isnan(sum(sum(gl)))==1 || isinf(sum(sum(gl)))==1
            disp(['t=', num2str(t)])
        end
        glt_1 = gl;
        Rt_1 = R;
        if PLM==1
            a = a + k(:,t).^(-1).*( ysim(:,t)-(a + b*xsim(:,t-1)) );
            gl = [a,b];
        elseif PLM==2
            R = R + k(:,t).^(-1)*([1;xsim(:,t-1)]*[1;xsim(:,t-1)]' - R);
            gl = (gl' + k(:,t).^(-1)*(R\([1;xsim(:,t-1)]*(ysim(:,t)-gl*[1;xsim(:,t-1)])')))';
        end
        % here I do a fake projection facility b/c can't take eig(gl) since
        % gl not square
        fake_eigs = eig((gl*gl').^(1/2));
        max_mod = max(abs(fake_eigs));
        if max_mod >= 1
            R = Rt_1;
            gl = glt_1;
        end
        fake_eigs_seq(:,:,t) = fake_eigs;
        max_mod_seq(t) = max_mod;
        gl_seq(:,:,t) = gl;
        R_seq(:,:,t) = R;
        
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
