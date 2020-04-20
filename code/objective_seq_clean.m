function [resids] = objective_seq_clean(seq,param,gx,hx,eta,PLM,gain,T,ndrop,e)
kapp = param.kapp;
sig  = param.sig;
psi_pi = param.psi_pi;
psi_x  = param.psi_x;

[s1, s2, s3, s4] = smat(param,hx);

% Simulate given input sequences
[xsim, ysim, ~, ~, FA, FB] = sim_learnLH_clean_given_seq(param,gx,hx,eta,PLM, gain, T,ndrop,e,seq);

% Evaluate residuals
pi = ysim(1,:);
x  = ysim(2,:);
i  = ysim(3,:);
s  = xsim;
fa = FA;
fb = FB;

resTR = -i + psi_pi*pi + psi_x*x;
resIS = -x -sig*i + s1*fb + s2*s;
resPC = -pi +kapp*x + s3*fa + s4*s;

resids = [resIS; resPC; resTR];


