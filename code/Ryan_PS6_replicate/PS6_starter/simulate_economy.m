function [longsim,LAM] = simulate_economy(Xgrid,longsim_init,pol_final,param,set)

load state_idx 
load v_idx
param_unpack


nt = size(longsim_init,2);
longsim = longsim_init;
H_final = pol_final;

for jj = 1:nt
    K = exp(longsim(state_idx(1),jj));
    GAM = exp(longsim(state_idx(2),jj));

    H = ndim_simplex_eval(Xgrid,longsim(state_idx,jj),H_final);
    
    C = (1-alph)/chi*GAM.^(alph/(alph-1)).*(K./H).^alph;
    K_p = GAM.^(alph/(alph-1)).*K.^alph.*H.^(1-alph) + (1-del)*K.*GAM.^(1/(alph-1)) - C;
    
    W = (1-alph)*GAM.^(alph/(alph-1)).*(K./H).^alph;
    R = alph*GAM*(K./H).^(alph-1);
    I =  K_p - (1-del)*K*GAM.^(1/(alph-1));

    %Assign values to longsim
    longsim(c_idx,jj) = C;
    longsim(h_idx,jj) = H;
    longsim(i_idx,jj) = I;
    longsim(r_idx,jj) = R;
    longsim(w_idx,jj) = W;
    if jj<nt
        longsim(state_idx(1),jj+1) = log(K_p);
    end
    
end

%These things compute outside of loop