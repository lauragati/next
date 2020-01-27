function EL = objective_CB(varp,setp,eN,burnin,PLM,gain)
% 26 Jan 2020

sig_r = setp.sig_r;
sig_i = setp.sig_i;
sig_u = setp.sig_u;

% an admittedly awkward structure to sub in the variable params
 setp.('psi_pi') = varp(1);
 setp.('psi_x') = varp(2);

param = setp;

% RE model
[fyn, fxn, fypn, fxpn] = model_NK(param);
[gx,hx]=gx_hx_alt(fyn,fxn,fypn,fxpn);
[ny, nx] = size(gx);
[Aa, Ab, As] = matrices_A_13_true_baseline(param, hx);
SIG = eye(nx).*[sig_r, sig_i, sig_u]';

[~,T,N] = size(eN);
y= zeros(ny,T,N);
for n=1:N
    e = squeeze(eN(:,:,n));
    [~, y(:,:,n)] = sim_learnLH(gx,hx,SIG,T,burnin,e, Aa, Ab, As, param, PLM, gain);
end

EL = loss_CB(param,y);
