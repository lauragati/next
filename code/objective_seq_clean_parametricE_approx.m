function [resids] = objective_seq_clean_parametricE_approx(seq,B,n_input_jumps,param,gx,hx,eta,PLM,gain,T,ndrop,e, alph_opt,x,fegrid, g_fe, knowTR)
% a version of objective_seq_clean with an expectation equation as a
% residual equation
% this version uses the approximated LOMgain from the GMM estimation, hence
% "approx"
% NOTE: this is for the univariate LOM gain, I only implement PEA for that.
% 13 June 2020
kapp = param.kapp;
sig  = param.sig;
% psi_pi = param.psi_pi;
% psi_x  = param.psi_x;
% rho_k = param.rho_k;
% gam_k = param.gam_k;
lamx = param.lamx;
alph= param.alph;
bet = param.bet;

mpshock=0;
[s1, s2, s3, s4] = smat(param,hx, knowTR,mpshock);

% Simulate given input sequences
% input k
% [xsim, ysim, k, ~, FA, FB, FEt_1] = sim_learnLH_clean_given_seq2(param,gx,hx,eta,PLM, gain, T,ndrop,e,seq,n_input_jumps);
% input FE
% [xsim, ysim, k, phi_seq, FA, FB, FEt_1,g_pi,g_pibar] = sim_learnLH_clean_given_seq3(param,gx,hx,eta,PLM, gain, T,ndrop,e,seq,n_input_jumps);
[xsim, ysim, k, phi_seq, FA, FB, FEt_1,g_pi] = sim_learnLH_clean_given_seq3_approx(param,gx,hx,eta,PLM, gain, T,ndrop,e,seq,n_input_jumps,...
    alph_opt,x,fegrid, g_fe, knowTR);

% Evaluate residuals (leave out first and last periods)
pi = ysim(1,2:end-1);
x  = ysim(2,2:end-1);
% xl = ysim(2,1:end-2);
i  = ysim(3,2:end-1);
s  = xsim(:,2:end-1);
st_1 = xsim(:,1:end-2);
fa = FA(:,2:end-1);
fb = FB(:,2:end-1);
fe = FEt_1(1,2:end-1); % select the FE(pi);
k1 = 1./k(2:end-1);
pibar = squeeze(phi_seq(1,1,1:end-2))';
b = gx*hx;
g_pi = g_pi(2:end-1);
% g_pibar = g_pibar(2:end-1);

% resTR = -i + psi_pi*pi + psi_x*x;
resIS = -x -sig*i + s1*fb + s2*s;
resPC = -pi +kapp*x + s3*fa + s4*s;
% resA6 = -k1(2:end) + rho_k*k1(1:end-1) + gam_k*fe.^2;
resA7 = -fe + pi - pibar - b(1,:)*st_1;
% resRETC = pi+lamx/kapp*x;
% res_commitTC = pi+lamx/kapp*x - lamx/kapp*xl;
% resanchTC0 = pi+lamx/kapp*x - lamx/kapp*(1-alph)*bet/(1-alph*bet)*(k1(2:end) + fe.*g_pi);
% Create basis as 1st, 2nd and 3rd powers of the states
X = [k1;pibar;s([1,end],:)];
sx = [X; X.^2; X.^3];
E = B'*sx;
resanchTC = pi+lamx/kapp*x - lamx/kapp*(1-alph)*bet/(1-alph*bet)*(k1 + fe.*g_pi).*E; 

resids = [resIS; resPC; resanchTC; resA7];

if min(k1)<0
    resids = ones(size(resids))*1e+10;
    disp('encountered negative k1, setting resids to 1e+10')
end








