% simulate data from learning model with Euler equation learning
% 12 sept 2019

function [z, conv] = sim_learn_EE(param,gx,R0,phi0,Ap,As,s,tol)

l = size(gx,1);
% Initialization
Rt = R0;
phit = phi0;
z = zeros(l,T);
h=1;
diff=zeros(T,1);
for t=1:T
    st = s(:,t);
    if t==1
        diff(t)=1;
        z(:,t)=gx*st;
    else
        % One-period-ahead expectation today
        Ezh = Ez_h(param, set, phit, st, h);
        
        % ALM
        z(:,t) = Ap*Ezh + As*st;
        
        % RLS
        Rp = Rt + t^(-1)*(st*st' - Rt_1);
        phip = phit + t^(-1)*Rp^(-1)*st*(z(:,t) - (phit+st));  
        
        % check convergence
        diff(t) = max(max(abs(phip - phit)))';
        if t>2 && diff(t) < tol && conv ==0
            conv = t;
        end
        
        Rt = Rp;
        phit = phip;

    end
end

