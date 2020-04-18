% computes the matrices I named "s" in front of fb, st, fa and st in
% equations A9 and A10 (e.g. materials25, appendix)
function [s1, s2, s3, s4, s5] = smat(param,hx)
kapp = param.kapp;
bet  = param.bet;
alph = param.alph;
sig  = param.sig;
nx = length(hx);

s1 = [sig,1-bet,-sig*bet];
s2 = sig*[1,0,0]*(eye(nx)-bet*hx)^(-1);
s3 = [(1-alph)*bet,kapp*alph*bet,0];
s4 = [0,0,1]*(eye(nx)-alph*bet*hx)^(-1);
% s5 = [0,1,0]; % if you wanna include a monpol shock.
s5 = [0,0,0]; % if you don't wanna include a monpol shock.


psi_pi = param.psi_pi;
psi_x  = param.psi_x;
s1_TR = [sig-sig*bet*psi_pi, 1-bet-sig*bet*psi_x,0];

% Impose condition (*): input s1_old for s1, which really is the
% assumption of agents knowing the Taylor rule:
s1 = s1_TR;
