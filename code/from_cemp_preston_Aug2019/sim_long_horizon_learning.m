function [z, conv, diff, a, b, FA, FB, FY] = sim_long_horizon_learning(param,R0,phi0,A1,A2,A3,s,tol)
% Evolution of economy, Preston
l = 3; % l is the number of regressands
[n,T] = size(s);
k =n+1;
w_seq = [ones(1,T);s]; %(k,T), these are the regressors

% Initialization
Rt_1 = R0;
phit_1 = phi0;

% add period t shock to nominal int rate
A4 = [0,0;0,0;1,0];

conv = 0;
z = zeros(l,T);
a = zeros(l,T);
b = zeros(l,n,T);
FA = zeros(k,T);
FB = zeros(k,T);
FY = zeros(n,T);
diff = zeros(T,1);
for t=1:T
    % Split phi(t-1) up into a(t-1) and b(t-1) neatly
    at_1 = phit_1(1,:)';
    bt_1 = phit_1(2:end,:)';
    
    % ALM
%     z(:,t) = ( (1-alph*bet)^(-1)*A1 + (1-bet)^(-1)*A2 )*at_1 ...
%         + ( A1*bt_1*(eye(n)-alph*bet*C)^(-1) + A2*bt_1*(eye(n)-bet*C)^(-1) + A3*f'*(eye(n)-bet*C)^(-1) ) *s(:,t) ...
%         + A4*s(:,t); % the last term is the effect of ibar_t on i

%     z(:,t) = ( 1/(1-alph*bet)*A1 + 1/(1-bet)*A2 )*at_1 ...
%         + ( A1*bt_1/(eye(n)-alph*bet*C) + A2*bt_1/(eye(n)-bet*C) + A3*f'/(eye(n)-bet*C) ) *s(:,t) ...
%         + A4*s(:,t); % the last term is the effect of ibar_t on i
    [fa, fb, fy] = fafbfy(at_1,bt_1,param,s(:,t));
    z(:,t) = A1*fa + A2*fb + A3*fy + A4*s(:,t);
    
    
    % Shock period t-1
    if t==1
        Rt = Rt_1;
        phit = phit_1;
    else
        wt_1 = w_seq(:,t-1);
        % RLS
        Rt = Rt_1 + t^(-1)*(wt_1*wt_1' - Rt_1);
%         phit = phit_1 + t^(-1)*Rt^(-1)*wt_1*(z(:,t-1) - phit_1'*wt_1)'; % z(t-1)  
        phit = phit_1 + t^(-1)*Rt^(-1)*wt_1*(z(:,t) - phit_1'*wt_1)'; %  z(t) <--
%         phit = phit_1 + t^(-1)*Rt^(-1)*wt_1*(z(:,t)' - wt_1'*phit_1); % this is identical to the previous
    end
    
    % check convergence
    at = phit(1,:)';
    bt = phit(2:end,:)';
    diff(t) = max(max(abs([at - at_1, bt - bt_1])))';
    if t>2 && diff(t) < tol && conv ==0
        conv = t;
    end
    
    % update
    Rt_1 = Rt;
    phit_1 = phit;
    % store
    a(:,t) = at;
    b(:,:,t) = bt;
    FA(:,t) = fa;
    FB(:,t) = fb;
    FY(:,t) = fy;
end

