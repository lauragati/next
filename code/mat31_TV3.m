function [tv, pibar, k1] = mat31_TV3(param,gx,hx,pp,i,pibart_1,k1t_1, s,sgrid,PI)
% Equations A9 and A10 in materials 25.
bet =param.bet;
lamx = param.lamx;
rhok = param.rho_k;
gamk = param.gam_k;
ns = length(sgrid);
b = gx*hx;

% use existing technology to compute these
[fa, fb] = fafb_anal_constant_free(param,[pibart_1;0;0],b,s,hx);

z = A9A10(param,hx,fa,fb,s,i);
pi = z(1);
x  = z(2);

% period loss
Lt =  pi^2 +lamx*x^2;

% Implied endogenous states at t
% crazy: compute forward and backward expectations at once
v = zeros(ns,ns,ns,ns);
% First: compute "expected fe" given what yesterday's shocks were
for i=1:ns
    rt_1=sgrid(i);
    for j=1:ns
        ut_1 = sgrid(j);
        st_1 = [rt_1;0;ut_1];
        fe = pi - (pibart_1+b(1,:)*st_1);
        k1 = rhok*k1t_1 + gamk*(fe)^2;
        pibar = pibart_1 + k1*(fe);
        
        % Second: compute expected future shocks for each expected past
        v = zeros(ns,ns);
        for ip=1:ns
            rp=sgrid(ip);
            for jp=1:ns
                up = sgrid(jp);
                v(i,j) = fnval(pp,{k1,pibar,rp,up});
            end
        end
        v(i,j,ip,jp) = fnval(pp,{k1,pibar,rt_1,ut_1});
    end
end

p = PI(1,1);

tv = Lt + bet*(v(:)'*p*ones(size(v(:))));

