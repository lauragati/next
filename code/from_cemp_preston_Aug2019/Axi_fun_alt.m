function Axi = Axi_fun_alt(gam, Gam, rho, fk)

Axi = [0,0,0;  0,rho,0; (1-gam)*Gam*fk^(-1),rho, gam ];