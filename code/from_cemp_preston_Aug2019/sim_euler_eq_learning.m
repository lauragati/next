function [z, conv, diff, a, b] = sim_euler_eq_learning(param,gx,R0,phi0,A1,A2,s,tol)
% Evolution of economy, Bullard & Mitra: (follows my notes on the paper)
f = param.f;
l = 3; % l is the number of regressands
[n,T] = size(s);
w_seq = [ones(1,T);s]; %(k,T), these are the regressors

% Initialization
Rt_1 = R0;
phit_1 = phi0;
z = zeros(l,T);

a = zeros(l,T);
b = zeros(l,n,T);

conv = 0;
diff = zeros(T,1);
for t=1:T
    if t==1
        diff(t)=1;
        z(:,t)=gx*s(:,t);
    else
        % Split phi(t-1) up into a(t-1) and b(t-1) neatly
        at_1 = phit_1(1,:)';
        bt_1 = phit_1(2:end,:)';
        
        % ALM
        z(:,t) = A1*at_1 + (A1*bt_1 + A2*f')*s(:,t);
        
        % RLS
        wt_1 = w_seq(:,t-1);
        Rt = Rt_1 + t^(-1)*(wt_1*wt_1' - Rt_1);
%         phit = phit_1 + t^(-1)*Rt^(-1)*wt_1*(z(:,t-1) - phit_1'*wt_1)';  % z(t-1) 
        phit = phit_1 + t^(-1)*Rt^(-1)*wt_1*(z(:,t) - phit_1'*wt_1)';  % z(t) <---
        
        % check convergence
        at = phit(1,:)';
        bt = phit(2:end,:)';
        diff(t) = max(max(abs([at - at_1, bt - bt_1])))';
        if t>2 && diff(t) < tol && conv ==0
            conv = t;
        end
        
        Rt_1 = Rt;
        phit_1 = phit;
        a(:,t) = at;
        b(:,:,t) = bt;
    end
end
