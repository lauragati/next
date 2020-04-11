function [resids] = objective_seq_fsolve(seq,param,gx,hx,eta,T,ndrop,e)
kapp = param.kapp;
sig  = param.sig;
psi_pi = param.psi_pi;
psi_x  = param.psi_x;

[stuff1, stuff2, stuff3, stuff4] = stuff(param,hx);

% Simulate given input sequences
[xsim, ysim, ~, ~, FA, FB] = sim_learnLH_clean_smooth_given_seq(param,gx,hx,eta,seq,T,ndrop,e);

% Evaluate residuals
pi = ysim(1,:);
x  = ysim(2,:);
i  = ysim(3,:);
s  = xsim;
fa = FA;
fb = FB;

resA3 = -i + psi_pi*pi + psi_x*x;
resA9 = -x -sig*i + stuff1*fb + stuff2*s;
resA10 = -pi +kapp*x + stuff3*fa + stuff4*s;

resids = [resA3; resA9; resA10];
% yo - I don't need to select which residuals I pick b/c if they are
% satisfied, so much the better!

