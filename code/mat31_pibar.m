function pibar = mat31_pibar(param,OM,pibart_1,k1,s,st_1,i)
Om5  = OM{5};
Om12 = OM{12};
kapp =param.kapp;
sig = param.sig;
b1 = OM{end};
pibar = pibart_1 + k1*(-kapp*sig*i +Om12*pibart_1 + Om5*s -b1*st_1);