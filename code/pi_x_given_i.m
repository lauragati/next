% Rewrite the general model w/o the Taylor-rule: given a sequence of i,
% what are x and pi? 
% See materials17.tex for a model summary and Notes 25 March 2020
% 25 March 2020
% Update 29 March 2020
function [pi_x] = pi_x_given_i(param,hx,fa,fb,s,iseq,t)
sig  = param.sig;
alph = param.alph;
bet  = param.bet;
kapp = param.kapp;
nx = size(hx,1);
% stuff1 = [sig,1-bet,-sig*bet];
stuff1 = [sig,1-bet,0];

stuff2 = sig*[1,0,0]*(eye(nx)-bet*hx)^(-1);
stuff3 = [(1-alph)*bet,kapp*alph*bet,0];
stuff4 = [0,0,1]*(eye(nx)-alph*bet*hx)^(-1);

capT = length(iseq);
dsi = 0;
for tau = t:capT
dsi = dsi + bet^(tau -t)*iseq(tau);
end

% pi_x = [0, 1; 1, -kapp]^(-1) * [-sig*iseq(t) + stuff1*fb + stuff2*s ; stuff3*fa + stuff4*s];

pi_x = [0, 1; 1, -kapp]^(-1) * [-sig*dsi + stuff1*fb + stuff2*s ; stuff3*fa + stuff4*s];
pi_x2 = [kapp, 1; 1, 0]      * [-sig*dsi + stuff1*fb + stuff2*s ; stuff3*fa + stuff4*s];

x = -sig*dsi + stuff1*fb + stuff2*s;
pit = kapp*x + stuff3*fa + stuff4*s;
% dbstop if warning
% if abs(sum(pi_x - [pit;x])) > 1e-4
%     pi_x - pi_x2
%     [pi_x(1) pit]
%     [pi_x(2), x]
%     warning('error')
% end

