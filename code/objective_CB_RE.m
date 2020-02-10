function EL = objective_CB_RE(varp,setp,eN,burnin)
% objective function of CB for RE NK model
% 9 Feb 2020

sig_r = setp.sig_r;
sig_i = setp.sig_i;
sig_u = setp.sig_u;

% an admittedly awkward structure to sub in the variable params
 setp.('psi_pi') = varp(1);
%  setp.('psi_x') = varp(2); % only optimize over psi_pi for now

param = setp;

% RE model
[fyn, fxn, fypn, fxpn] = model_NK(param);
[gx,hx]=gx_hx_alt(fyn,fxn,fypn,fxpn);
[ny, nx] = size(gx);
SIG = eye(nx).*[sig_r, sig_i, sig_u]';

[~,T,N] = size(eN);
y= zeros(ny,T,N);
for n=1:N
    e = squeeze(eN(:,:,n));
    [~, y(:,:,n)] = sim_model(gx,hx,SIG,T,burnin,e);
end

EL = loss_CB(param,y);
