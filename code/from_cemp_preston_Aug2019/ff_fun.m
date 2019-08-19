function f = ff_fun(fk, pibar_t, gam, Gam)

fpibar = fpibar_fun(fk,pibar_t,gam,Gam);
fxi = fxi_fun(gam, Gam, fpibar);

f = [fpibar, fxi']';