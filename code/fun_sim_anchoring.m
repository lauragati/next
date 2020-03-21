%fun_sim_anchoring.m
% Generate data from anchoring model
% A "function" copy of command_sim_anchoring.m
% 20 March 2020
function y = fun_sim_anchoring(param,T,N, burnin,eN,PLM,gain)
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

y = zeros(ny,T,N);
for n=1:N
    % Sequence of innovations
    e = squeeze(eN(:,:,n));
    
    % Learning
    [~, y_LH] = sim_learnLH(gx,hx,SIG,T+burnin,burnin,e, Aa, Ab, As, param, PLM, gain);
    y(:,:,n) = y_LH;
end