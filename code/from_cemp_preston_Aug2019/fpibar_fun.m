function fpibar = fpibar_fun(fk,pibar,gam,Gam)

fpibar = (1-(1-gam)*(1-Gam)*fk.^(-1)).*pibar;