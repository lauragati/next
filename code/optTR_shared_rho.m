function [psi_pi_opt, phi_pi_opt] =  optTR_shared_rho(rhoall,lami)
% Optimal TR coefficients on inflation for the case when all shocks are
% driven by the same AR-coefficient, rhoall.
% phi is RE. psi is learning
% 3 Feb 2020

param = parameters_next;
bet  = param.bet;
alph = param.alph;
sig  = param.sig;
kapp = param.kapp;

phi_pi_opt = kapp.*sig.*(lami.*((-1)+rhoall).*((-1)+bet.*rhoall)+(-1).*kapp.* ...
  lami.*rhoall.*sig).^(-1);

psi_pi_opt = kapp.*((-1)+bet.*((-1)+rhoall)).*((-1)+alph.*bet.*((-1)+rhoall)).* ...
  sig.*(bet.*lami.*((-1)+bet+alph.*bet.*((-1)+rhoall)).*((-1)+ ...
  rhoall)+kapp.*lami.*((-1)+alph.*bet.*((-1)+rhoall)).*sig).^(-1);