% MODEL_SS - Return the steady state of the model computed analytically
%
% usage:
% 
% [ss, parameters] =model_ss(param)


function [ss,param,set] = model_ss(param,set)

set.empty = [];

%Upack parameters object
param_unpack

%BEGIN_EXTRACT_HERE

%Use closed form expressions for the ss values.
r = (1/(bet*gam^(-1/(1-alph))))-1+del;
kh = (r/(alph*gam))^(1/(alph-1));
w = (1-alph)*gam^(alph/(alph-1))*(kh)^alph;
ik = (gam^(1/(1-alph))-(1-del))/(gam^(1/(1-alph)));
c = w/chi;
ck = gam^(alph/(alph-1))*(kh)^(alph-1)-ik;

k = ck^-1*c;
i = ik*k;
h = kh^-1*k;

%Put the ss values in a vector consistent with Y and X vectors in model.m
Yss  = [c h w r i];
Xss  = [k gam];
ss  = [Yss Xss];


%END_EXTRACT_HERE


set.kbar = k;
set.cbar = c;
set.rbar = r;
set.hbar = h;