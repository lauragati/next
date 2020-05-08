function [resids] = objective_seq_clean2(seq,param,gx,hx,eta,PLM,gain,T,ndrop,e)
% a version of objective_seq_clean with an expectation equation as a
% residual equation
kapp = param.kapp;
sig  = param.sig;
psi_pi = param.psi_pi;
psi_x  = param.psi_x;
rho_k = param.rho_k;
gam_k = param.gam_k;

[s1, s2, s3, s4] = smat(param,hx);

% Simulate given input sequences
[xsim, ysim, k, ~, FA, FB, FEt_1] = sim_learnLH_clean_given_seq(param,gx,hx,eta,PLM, gain, T,ndrop,e,seq);

% Evaluate residuals (leave out first and last periods)
pi = ysim(1,2:end-1);
x  = ysim(2,2:end-1);
i  = ysim(3,2:end-1);
s  = xsim(:,2:end-1);
fa = FA(:,2:end-1);
fb = FB(:,2:end-1);
fe = FEt_1(1,2:end-1); % select the FE(pi);
k1 = 1./k(1:end-1);

resTR = -i + psi_pi*pi + psi_x*x;
resIS = -x -sig*i + s1*fb + s2*s;
resPC = -pi +kapp*x + s3*fa + s4*s;
res_gain = -k1(2:end) + rho_k*k1(1:end-1) + gam_k*fe.^2;

% resids = [resIS; resPC; resTR];
resids = [resIS; resPC; resTR; res_gain];
% resids=resTR;


