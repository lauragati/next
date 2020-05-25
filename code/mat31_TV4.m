function Ln = mat31_TV4(param,gx,hx,i,pibart_1,s,v)
% once you have the value function for a sequence of shocks, find optimal i-sequence
bet =param.bet;
lamx = param.lamx;
sig  = param.sig;
kapp = param.kapp;
b = gx*hx;

a = [pibart_1;zeros(size(pibart_1));zeros(size(pibart_1))];
% use existing technology to compute these
[fa, fb] = fafb_anal_constant_free(param,a,b,s,hx);


% implement A9A10 for matrices
[s1, s2, s3, s4] = smat(param,hx);
x = -sig.*i + s1*fb + s2*s;
pi = kapp*x + s3*fa + s4*s;

% loss
L =  pi.^2 +lamx.*x.^2 + bet.*v;

% norm
Ln = max(max(abs(L)));


