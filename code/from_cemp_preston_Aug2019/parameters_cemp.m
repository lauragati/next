function [param0, pmean, pstd, pshape, set, variable, param_correct] = parameters_cemp
% pistar, thetbar, gbar, gamma, Gamma, rhophi, sige, sigmu, sigo1, sigo2, sigo3, sigo4,sigo5

% 1= beta; 2=gamma; 3 = normal; 4 = inv wishart; 5 = uniform; 6=inv gamma.
pshape = [3,2,2,1,1,1,6,6,6,6,6,6,6 ];
%Normal Gamma Gamma Beta Beta Beta Inv.-Gamma Inv.-Gamma Inv.-Gamma Inv.-Gamma Inv.-Gamma Inv.-Gamma Inv.-Gamma
pmean = [2.250 0.050 0.100 0.500 0.500 0.500 0.100 0.100 0.100 0.100 0.100 0.100 0.100];
pstd  = [0.400 0.040 0.090 0.265 0.265 0.200 2.000 1.000 1.000 1.000 1.000 1.000 1.000];
nparam= length(pshape);

% if nothing else specified, set param0 to their prior mean
param0 = pmean;

param_correct =  [2.472 0.029     0.145 0.128  0.891  0.877   0.084 0.359  0.277 0.042   0.021  0.073 0.049];
set.pistar  = param_correct(1);
set.thetbar = param_correct(2);
set.gbar    = param_correct(3);
set.gam     = param_correct(4);
set.Gam     = param_correct(5);
set.rho     = param_correct(6);
set.sige    = param_correct(7);
set.sigmu   = param_correct(8);
set.o1      = param_correct(9);
set.o2      = param_correct(10);
set.o3      = param_correct(11);
set.o4      = param_correct(12);
set.o5      = param_correct(13);

% if variable(3)=1, that means gbar belongs to param, if it's 0, then gbar
% belongs to set
% variable = [1,1,0,1,1,1,1,1,  1,1,1,1,1];
variable = [1,1,1,1,1,1,1,1,  1,1,1,1,1];
% variable = [0,0,0,0,0,0,0,0,  0,0,0,0,0];