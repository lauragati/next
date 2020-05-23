function k1 = mat31_k1(OM,param,k1t_1,i,pibart_1,s,st_1)
Om5  = OM{5};
Om12 = OM{12};
kapp =param.kapp;
sig = param.sig;
rhok = param.rho_k;
gamk = param.gam_k;
b1 = OM{end};

k1 = rhok*k1t_1 + gamk*(-kapp*sig*i +Om12*pibart_1 + Om5*s -b1*st_1)^2;