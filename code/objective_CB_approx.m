function EL = objective_CB_approx(varp,setp,eN,ndrop,PLM,gain,alph,x, knowTR)
% 30 June 2020

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

% knowTR=1;
mpshock=1;

[~,T,N] = size(eN);
y= zeros(ny,T,N);
for n=1:N
    e = squeeze(eN(:,:,n));
    [~, y(:,:,n)] = sim_learnLH_clean_approx(alph,x,param,gx,hx,SIG, PLM, gain, T,ndrop,e,knowTR,mpshock);
    % should use ... univariate.m to use univariate anchoring function, but there may be a problem with this code...
%     [~, y(:,:,n)] = sim_learnLH_clean_approx_univariate(alph,x,param,gx,hx,SIG, PLM, gain, T,ndrop,e,v,knowTR,mpshock);
end

EL = loss_CB(param,y);
