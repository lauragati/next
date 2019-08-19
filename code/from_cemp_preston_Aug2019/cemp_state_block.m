function [fk, fpibar, Apibar, fxi, Axi] = cemp_state_block(param, k_t_1, pibar_t_1)
% 22 June 2019
%
% All of this given particle i

pistar  = param(1);
thetbar = param(2);
gbar    = param(3);
gam     = param(4);
Gam     = param(5);
rho     = param(6);
sige    = param(7);
sigmu = param(8);
sigo1 = param(9);
sigo2 = param(10);
sigo3 = param(11);
sigo4 = param(12);
sigo5 = param(13);

% State block
% fk
I = abs((1-gam)*(Gam-1)*pibar_t_1) <= thetbar*(sige+sigmu);
% disp(['Theta is: ', num2str(abs((1-gam)*(Gam-1)*pibar_t_1)), '; threshold is: ', num2str(thetbar*(sige+sigmu))])

fk = I.*(k_t_1+1) + (1-I).*gbar^(-1);

% fpibar
fpibar = (1-(1-gam)*(1-Gam)*fk.^(-1)).*pibar_t_1;

% Apibar
Apibar = [fk^(-1), 0, 0];

% fxi
fxi = [0 ; 0; (1-gam)*Gam*fpibar];

% Axi
Axi = vertcat(zeros(1,3), [0 rho, 0], [(1-gam)*Gam*fk^(-1),rho, gam]);