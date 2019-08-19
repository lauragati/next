function A = A_fun(fk, gam, Gam, rho)
Apibar = Apibar_fun(fk);
Axi = vertcat(zeros(1,3), [0 rho, 0], [(1-gam)*Gam*fk^(-1),rho, gam]);
A = [Apibar; Axi];