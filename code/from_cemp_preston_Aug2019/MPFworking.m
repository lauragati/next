function [LL, xn_swarm, xl_swarm, Ptt_seq, qt_seq, qtilde_seq, K] = MPFworking(Y,param, N,obs_blocksize)
% Cleaned-up version of MPF_test_resampling1_1by1
% 8 August 2019

% Y is the whole sample, w/o missing values
T = size(Y,2);
% l = size(Y,1);

k_n = 2; % the number of nonlinear states
k_l = 3; % the number of linear states
% ne = 2; % number of state innovations

[~, ~, ~, ~, set, variable] = parameters_cemp;
% Here I choose which parameters are estimated (which are variable, not
% constant).
% The idea is to set everything to the constant values 'set', and then
% replace the variable ones with the input parameter values.

variables = find(variable);
if isempty(variables)==1
    disp('not estimating any params')
    paramset = param;
else
    paramset = struct2array(set);
    paramset(variables) = param;
end

gbar    = paramset(3);
gam     = paramset(4);
rho     = paramset(6);
sige    = paramset(7);
sigmu = paramset(8);

% Things that are data- and particle-independent
[h0, hpibar, H, R] = cemp_obs_block(paramset,obs_blocksize);
[SIG, Sxi] = cemp_SIG_S(param); 


% 1.) Initialization
randgbar = gbar.*rand(1,N);
k10 = randgbar.^(-1);
pibar10 = randn(1,N);
xi10 = repmat([0, 0, 0]',1,N);
P10  = eye(k_l).*[sige+sigmu, sige/(1-rho^2), (sige+sigmu)/(1-gam^2)]; % this is a question! but not essential for dynamics

w10 = 1/N*ones(N,1); % initial resampled-weights

% Initialize output needed for smoother
qt_seq = zeros(T,N);
qtilde_seq = zeros(T,N);
xn_swarm = zeros(T,k_n,N);
xl_swarm = zeros(T,k_l,N);
Ptt_seq  = zeros(T, k_l, k_l);

% Initialize filter-internal things
k_t_1 = k10;
pibar_t_1 = pibar10;
xi_t_1 = xi10;
P_t_1 = P10;
w_t_1 = w10;

xip_i = zeros(k_l,1);
xi_t = zeros(k_l,N);
kp = zeros(1,N);
pibarp = zeros(1,N);
LL = -T*log(N);

for t=1:T
    % 2.) Importance weights
    Om_t = H'*P_t_1*H + R;
    MU = h0 + hpibar*pibar_t_1 + H'*xi_t_1;
    SIGM = Om_t;
        qt = w_t_1.* mvnpdf(Y(:,t)',MU',SIGM); % (N x 1)
    sumqt = sum(qt);
    qtilde_t = qt./sumqt;
    sumqtilde2 = sum(qtilde_t.^2);
    ESS = 1/sumqtilde2;
    % 3.) Resampling
    if sumqt ~= 0 && ESS < 0.75*N
        indeces = randsample(1:N,N,'true',qtilde_t);
        pibar_t = pibar_t_1(indeces);
        k_t  = k_t_1(indeces);
        w_t =1/N*ones(N,1);
    else % don't resample
        pibar_t = pibar_t_1;
        k_t  = k_t_1;
        w_t  = qtilde_t;
    end
    % 4) Linear measurement equation / Kalman filter inference
    K = P_t_1*H*Om_t^(-1);
    P_t = P_t_1-K*H'*P_t_1;
%     P_t = P_t_1-K*Om_t*K'; % this and the previous line are identical
    Peta = P_t(1,1); Pphi = P_t(2,2);Ppi = P_t(3,3); Petaphi = P_t(1,2); Petapi = P_t(1,3); Pphipi = P_t(2,3);
    for i=1:N
        xi_t(:,i) = xi_t_1(:,i) + K*(Y(:,t) - h0 -hpibar*pibar_t(i) - H'*xi_t_1(:,i));
        % 5.) Particle filter prediction
        [fk_i, fpibar_i, Apibar_i, fxi_i, Axi_i] = cemp_state_block(paramset, k_t(i), pibar_t(i));
        kp(i) = fk_i;
        z_i = xi_t(1,i);
        MU_5 = fpibar_i + fk_i.^(-1).*z_i; 
        SIGM_5 = fk_i.^(-2)*Peta;
        pibarp(i) = mvnrnd(MU_5,SIGM_5,1); %  3.) let matlab scale (that is scale with chol(SIG))

        % 6.) Linear model prediction
        ktilde_i = P_t*Apibar_i'*(Apibar_i*P_t*Apibar_i')^(-1);
        xitilde_i = xi_t(:,i) + ktilde_i*(pibarp(i)-fpibar_i-fk_i.^(-1)*z_i);
        xip_i(:,i) = fxi_i + Axi_i*xitilde_i;
    end
    xip=xip_i;
    i1 = (Pphi - (Petaphi^2)/(Peta))*rho^2;
    i2 = (-(Petaphi*Petapi)/(Peta) - Pphipi)*rho*gam;
    Ptilde = [0, 0,     0;...
        0, i1,    i1+i2;...
        0, i1+i2, 2*i2+i1+(Ppi - (Petapi^2)/(Peta))*gam^2 ];
    
    Pp = Sxi*SIG*Sxi' + Ptilde;
    
    % update
    pibar_t_1 = pibarp;
    k_t_1  = kp;
    xi_t_1 = xip;
    w_t_1 = w_t;
    P_t_1 = Pp;
    
    % Store for smoother
    qt_seq(t,:) = qt';
    qtilde_seq(t,:) = qtilde_t';
    xn_swarm(t,:,:) = [k_t; pibar_t];
    xl_swarm(t,:,:) = xi_t;
    Ptt_seq(t,:,:) = P_t;
    
    %avoid adding log(0) -> only add it if sumqt~=0
    if sumqt ~=0
        LL = LL + log(sumqt);
    end
end