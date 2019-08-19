function Axi = Axi_fun(gam, Gam, rho, fk,m)

top_row = zeros(1,3,m);
middle_row = repmat([0 rho 0],1,1,m);
bottom_row = [(1-gam)*Gam*fk.^(-1)', rho*ones(size(fk))', gam*ones(size(fk))']';
Axi = cat(1,top_row, middle_row);
Axi(3,:,:) = bottom_row;