function [tv, pibar, k1] = mat31_TV3(param,gx,hx,pp,i,pibart_1,k1t_1,s,st_1,sgrid,PI)
% Equations A9 and A10 in materials 25.
bet =param.bet;
lamx = param.lamx;
rhok = param.rho_k;
gamk = param.gam_k;
rho_r = param.rho_r;
rho_u = param.rho_u;
sig_r = param.sig_r;
sig_u = param.sig_u;
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
% First: compute "expected fe" given what yesterday's shocks were
fe = pi - (pibart_1+b(1,:)*st_1);
k1 = rhok*k1t_1 + gamk*(fe)^2;
pibar = pibart_1 + k1*(fe);
% Second: compute expected future shocks for each current state
v = zeros(ns,ns);
r = s(1);
u = s(3);
for i=1:ns
    rp=rho_r*r + sig_r*sgrid(i);
    for j=1:ns
        up = rho_u*u +sig_u*sgrid(j);
        v(i,j) = fnval(pp,{k1,pibar,r,u,rp,up});
    end
end


tv = Lt + bet*(v(:)'*PI(:));
