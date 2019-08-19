function param = param_preston

n = 2;
f = eye(n,n);
rho = 0.9;
C = rho*eye(n);
alph = 0.059; % 0.6, or 0.059 to give a kappa of approx 15
bet  = 0.98;
sig  = 1; %0.5;
% kapp = 15; % 15
kapp = (1-alph)*(1-alph*bet)/alph;
psi_x  = 1.5; %0.9, 1.5
psi_pi = 1.5; %1.5
w_small = 1 +sig*psi_x +kapp*sig*psi_pi;

% Check some conditions

% Taylor principle (Preston, p. 110 (Mac 31). )
taylor_principle = kapp*(psi_pi-1)+(1-bet)*psi_x > 0;

eig(C) < 1;

param.n = 2;
param.f = f;
param.C = C;
param.alph = alph;
param.bet = bet;
param.sig = sig;
param.kapp = kapp;
param.psi_x = psi_x;
param.psi_pi = psi_pi;
param.w_small = w_small;