%fun_sim_anchoring_given_i.m
% Generate data from anchoring model, given an exog i-sequence
% A near-copy of fun_sim_anchoring.m
% 25 March 2020
function [pi,x,k,pibar,b,s,g_pi] = fun_sim_anchoring_given_i(param,capT,N,burnin,eN,PLM,gain,i_seq)
sig_r = param.sig_r;
sig_i = param.sig_i;
sig_u = param.sig_u;

% RE model
[fyn, fxn, fypn, fxpn] = model_NK(param);
[gx,hx]=gx_hx_alt(fyn,fxn,fypn,fxpn);
[ny, nx] = size(gx);
[Aa, Ab, As] = matrices_A_13_true_baseline(param, hx);
SIG = eye(nx).*[sig_r, sig_i, sig_u]';
eta = SIG; %just so you know

%% Simulate model

pi = zeros(capT,N);
x = zeros(capT,N);
k = zeros(capT,N);
pibar = zeros(capT,N);
b = zeros(ny,nx);
s = zeros(nx,capT,N);
g_pi = zeros(capT,N);
for n=1:N
    % Sequence of innovations
    e = squeeze(eN(:,:,n));
    
    % Learning
    [pi(:,n),x(:,n),k(:,n),pibar(:,n),b,s(:,:,n),g_pi(:,n)] = sim_learnLH_given_i(gx,hx,SIG,capT+burnin,burnin,e, Aa, Ab, As, param, PLM, gain, i_seq);
end