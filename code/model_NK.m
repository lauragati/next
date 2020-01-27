% 3-equation NK model with 3 shocks: natural rate, monetary policy and
% cost-push
% 12 September 2019
% adapted from model_preston.m

function [fyn, fxn, fypn, fxpn] = model_NK(param)

%Steady State
[ss, param] = model_NK_ss(param);

%Declare parameters
bet = param.bet;  
sig = param.sig;
alph = param.alph;
kapp = param.kapp;
psi_x = param.psi_x;
psi_pi = param.psi_pi;
w = param.w;
gbar = param.gbar;
thetbar = param.thetbar;
rho_r = param.rho_r;
rho_i = param.rho_i;
rho_u = param.rho_u;
sig_r = param.sig_r;
sig_i = param.sig_i;
sig_u = param.sig_u;

%Declare Needed Symbols
syms RN RN_p IB IB_p U_p U
syms PI PI_p XX XX_p I I_p

%Declare X and Y vectors
X  = [RN IB U]; 
XP = [RN_p IB_p U_p];

Y  = [PI XX I  ];
YP = [PI_p XX_p I_p ] ;


%Model Equations (materials2 and 3)
f(1) = -XX + XX_p - sig*(I - PI_p)+ sig*RN;
f(2) = -PI + kapp*XX +bet*PI_p +U; 
f(3) = -I + psi_pi*PI + psi_x*XX + IB;
f(4) = IB_p - rho_i*IB;
f(5) = RN_p - rho_r*RN;
f(6) = U_p - rho_u*U;

%Check Computation of Steady-State Numerically
fnum = double(subs(f, [Y X YP XP], [ss, ss]));
% disp('Checking steady-state equations:')
% disp(fnum);

%Log-linear approx
log_var = [];
f = subs(f, log_var, exp(log_var)); 
   
%Differentiate
fx  = jacobian(f,X);
fy  = jacobian(f,Y);
fxp = jacobian(f,XP);
fyp = jacobian(f,YP);

%Plug back into levels
f =   subs(f ,  log_var, log(log_var));
fx =  subs(fx , log_var, log(log_var));
fy =  subs(fy , log_var, log(log_var));
fxp = subs(fxp, log_var, log(log_var));
fyp = subs(fyp, log_var, log(log_var));

%Compute numerical values
fxn =  double(subs(fx , [Y X YP XP], [ss, ss]));
fyn =  double(subs(fy , [Y X YP XP], [ss, ss]));
fxpn = double(subs(fxp, [Y X YP XP], [ss, ss]));
fypn = double(subs(fyp, [Y X YP XP], [ss, ss]));
