function xn = F_fun(xn_1, vn)

param0  = parameters_cemp;
pistar  = param0(1);
thetbar = param0(2);
gbar    = param0(3);
gam     = param0(4);
Gam     = param0(5);
rhophi  = param0(6);
sige    = param0(7);
sigmu = param0(8);
sigo1 = param0(9);
sigo2 = param0(10);
sigo3 = param0(11);
sigo4 = param0(12);
sigo5 = param0(13);

Sxi = ones(3,2); % shock impact in lin state equation
Sxi(2,2) = 0;

k_t_1 = xn_1(1);
pibar_t_1 = xn_1(2);
xi_t_1 = xn_1(3:end);

fk = fk_fun(k_t_1,pibar_t_1,gam,Gam,sige,gbar,thetbar);
fpibar = fpibar_fun(fk,pibar_t_1,gam,Gam);
Apibar = Apibar_fun(fk);
Axi = Axi_fun(gam, Gam, rhophi, fk);
fxi = fxi_fun(gam, Gam, fpibar);

f = [fk; fpibar; fxi];
A = [zeros(1,3); Axi; Apibar];
S = [zeros(2,2); Sxi];
xn = f + A*xi_t_1 + S*vn;