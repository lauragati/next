function EL = objective_CB_approx(varp,setp,eN,ndrop,PLM,gain,alph,x, knowTR)
% 30 June 2020

sig_r = setp.sig_r;
sig_i = setp.sig_i;
sig_u = setp.sig_u;

% an admittedly awkward structure to sub in the variable params
setp.psi_pi = varp(1);
%  setp.('psi_x') = varp(2); % only optimize over psi_pi for now

param = setp;

% RE model
% try
[fyn, fxn, fypn, fxpn] = model_NK(param);
% catch
% end
[gx,hx]=gx_hx_alt(fyn,fxn,fypn,fxpn);
[ny, nx] = size(gx);
SIG = eye(nx).*[sig_r, sig_i, sig_u]';

% knowTR=1;
mpshock=1;

[~,T,N] = size(eN);
y= zeros(ny,T,N);
for n=1:N
    e = squeeze(eN(:,:,n));
    
    %     [~, y(:,:,n)] = sim_learnLH_clean_approx(alph,x,param,gx,hx,SIG, PLM, gain, T,ndrop,e,knowTR,mpshock);
    % should use ... univariate.m to use univariate anchoring function, but there may be a problem with this code...
    
    v = zeros(4,T);
    try
        [~, y(:,:,n), ~, ~, ~, ~, ~,~, explode_count, negk_count, explode_t, negk_t] = sim_learnLH_clean_approx_univariate(alph,x,param,gx,hx,SIG, PLM, gain, T,ndrop,e,v,knowTR,mpshock);
        
%         % Annualize inflation and interest rate - might not need to since it just scales up loss
%         y([1,3],:) = ((y([1,3],:)/100+1).^4 -1)*100;
       
        
    catch err
        disp(['Exploded (fe input to ndim_simplex was nan: history n = ', num2str(n)])
        fprintf(1,'The identifier was:\n%s',err.identifier);
        fprintf(1,'\n The error message was:\n%s',err.message);
        fprintf(1,'\n');
        y(:,:,n) = inf(ny,T);
    end
end

EL = loss_CB(param,y);
