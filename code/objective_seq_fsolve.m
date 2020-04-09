function [resids] = objective_seq_fsolve(seq,param,e,T,burnin,PLM,gain,gx,hx,SIG,Aa,Ab,As, H, variant)

% Evaluate residuals from simulation given input sequence(s)
[~,~,~,~,~,~,~,~,~,~,~,~,~,resids] = sim_learnLH_given_seq(gx,hx,SIG,T+burnin,burnin,e, Aa, Ab, As, param, PLM, gain, seq, H, variant);

