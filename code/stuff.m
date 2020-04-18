% computes the matrices I named "stuff" in front of fb, st, fa and st in
% equations A9 and A10 (e.g. materials25, appendix)
function [stuff1, stuff2, stuff3, stuff4, stuff5] = stuff(param,hx)
kapp = param.kapp;
bet  = param.bet;
alph = param.alph;
sig  = param.sig;
nx = length(hx);

stuff1 = [sig,1-bet,-sig*bet];
stuff2 = sig*[1,0,0]*(eye(nx)-bet*hx)^(-1);
stuff3 = [(1-alph)*bet,kapp*alph*bet,0];
stuff4 = [0,0,1]*(eye(nx)-alph*bet*hx)^(-1);
stuff5 = [0,1,0]; % if you wanna include a monpol shock.
% stuff5 = [0,0,0]; % if you don't wanna include a monpol shock.


% The below part verifies that it's stuff1 that is different between A9A10
% and the old A-matrices, because stuff1_old incorporates the assumption
% (*), which says that agents know the Taylor rule.
psi_pi = param.psi_pi;
psi_x  = param.psi_x;
w = param.w;
stuff1_old = [sig-sig*bet*psi_pi, 1-bet-sig*bet*psi_x,0];
Ab_old = [kapp/w; 1/w; (psi_x+kapp*psi_pi)/w]*stuff1_old;

% Impose condition (*): input stuff1_old for stuff1, which really is the
% assumption of agents knowing the Taylor rule:
stuff1 = stuff1_old;
