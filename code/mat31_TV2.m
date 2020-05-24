function [tv, pibar, k1] = mat31_TV2(param,gx,hx,pp,i,pibart_1,k1t_1, s, st_1,sgrid,PI)
% Equations 1-8 in materials 25.
alph =param.alph;
bet =param.bet;
kapp =param.kapp;
sig = param.sig;
lamx = param.lamx;
rhok = param.rho_k;
gamk = param.gam_k;
ns = length(sgrid);
nx = size(hx,2);
e1 = [1,0,0]; 
e3 = [0,0,1];
b = gx*hx;
b1 = b(1,:);
b2 = b(2,:);
b3 = b(3,:);

fa = 1/(1-alph*bet)*pibart_1 + b1*(eye(nx)-alph*bet*hx)^(-1)*s;
fb = 1/(1-bet)*pibart_1 + b1*(eye(nx)-bet*hx)^(-1)*s;
x = -sig*i + sig*fb + (1-bet)*b2*(eye(nx)-bet*hx)^(-1)*s -sig*bet*b3*(eye(nx)-bet*hx)^(-1)*s + sig*e1*(eye(nx)-bet*hx)^(-1)*s;
pi = kapp*x +(1-alph)*bet*fa +kapp*alph*bet*b2*(eye(nx)-alph*bet*hx)^(-1)*s +e3*(eye(nx)-alph*bet*hx)^(-1)*s;

% period loss
Lt =  pi^2 +lamx*x^2;

% Implied endogenous states at t
fe = pi - (pibart_1+b1*st_1);
k1 = rhok*k1t_1 + gamk*(fe)^2;
pibar = pibart_1 + k1*(fe);

% approximation of value function % here form weighted average
% expectations over tomorrow's r and u
v = zeros(ns,ns);
for i=1:ns
    rp=sgrid(i);
    for j=1:ns
        up = sgrid(j);
        v(i,j) = fnval(pp,{k1,pibar,rp,up});
    end
end

tv = Lt + bet*(v(:)'*PI(:));

