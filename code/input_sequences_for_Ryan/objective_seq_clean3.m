function [resids] = objective_seq_clean3(seq,n_input_jumps,param,gx,hx,eta,PLM,gain,T,ndrop,e)
% a version of objective_seq_clean with an expectation equation as a
% residual equation
kapp = param.kapp;
sig  = param.sig;
psi_pi = param.psi_pi;
psi_x  = param.psi_x;
rho_k = param.rho_k;
gam_k = param.gam_k;
lamx = param.lamx;
alph= param.alph;
bet = param.bet;

[s1, s2, s3, s4] = smat(param,hx);

% Simulate given input sequences
% input k
% [xsim, ysim, k, ~, FA, FB, FEt_1] = sim_learnLH_clean_given_seq2(param,gx,hx,eta,PLM, gain, T,ndrop,e,seq,n_input_jumps);
% input FE
[xsim, ysim, k, phi_seq, FA, FB, FEt_1,g_pi,g_pibar] = sim_learnLH_clean_given_seq3(param,gx,hx,eta,PLM, gain, T,ndrop,e,seq,n_input_jumps);

% Evaluate residuals (leave out first and last periods)
pi = ysim(1,2:end-1);
x  = ysim(2,2:end-1);
i  = ysim(3,2:end-1);
s  = xsim(:,2:end-1);
st_1 = xsim(:,1:end-2);
fa = FA(:,2:end-1);
fb = FB(:,2:end-1);
fe = FEt_1(1,2:end-1); % select the FE(pi);
k1 = 1./k(1:end-1);
pibar = squeeze(phi_seq(1,1,1:end-2))';
b = gx*hx;
g_pi = g_pi(2:end-1);
g_pibar = g_pibar(2:end-1);

resTR = -i + psi_pi*pi + psi_x*x;
resIS = -x -sig*i + s1*fb + s2*s;
resPC = -pi +kapp*x + s3*fa + s4*s;
resA6 = -k1(2:end) + rho_k*k1(1:end-1) + gam_k*fe.^2;
resA7 = -fe + pi - pibar - b(1,:)*st_1;
resRETC = pi+lamx/kapp*x;
resanchTC0 = pi+lamx/kapp*x - (1-alph)*bet/(1-alph*bet)*(k1(2:end) + fe.*g_pi);
% interesting: this first part of the anchTC is just a mu below resRETC.

% resids=resTR;
% resids = [resIS; resPC; resTR];
% resids = [resIS; resPC; resTR; resA6];
% resids = [resTR; resA6];
% resids = [resIS; resPC; resTR; resA7];
% resids = [resTR; resA7];
% resids = [resRETC; resA7];
% resids = [resanchTC0; resA7];
% resids = [resanchTC0];
resids = [resIS; resPC; resanchTC0; resA7];
% resids = [resIS; resanchTC0; resA7];







