function loss = objective_seq(seq,param,e,T,burnin,PLM,gain,gx,hx,SIG,Aa,Ab,As,y,s_r)

% Evaluate simulation given input sequence
[~, yg] = sim_learnLH_given_seq(gx,hx,SIG,T+burnin,burnin,e, Aa, Ab, As, param, PLM, gain, seq);

% Evaluate errors compared to Taylor rule simulation
resids = s_r.*(y-yg);

loss = max(max(abs(resids)));