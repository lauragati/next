% 3-equation NK model with 3 shocks: natural rate, monetary policy and
% cost-push
% 21 Jan 2020
% adapted from model_NK.m

function [fyn, fxn, fypn, fxpn] = model_NK_Markov_switchingTR(param)

%Steady State
[ss, param] = model_NK_Markov_switchingTR_ss(param);

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
psi1 = param.psi1;
psi2 = param.psi2;
p11 = param.p11;
p21 = param.p21;
p22 = param.p22;
p12 = param.p12;

%Declare Needed Symbols
syms RN RN_p IB IB_p U_p U
syms PI1 PI1_p XX1 XX1_p I1 I1_p
syms PI2 PI2_p XX2 XX2_p I2 I2_p


%Declare X and Y vectors
X  = [RN IB U]; 
XP = [RN_p IB_p U_p];

Y  = [PI1 PI2 XX1 XX2 I1 I2  ];
YP = [PI1_p PI2_p XX1_p XX2_p I1_p I2_p] ;


%Model Equations (materials2 and 3)
f(1)     = -XX1 + (p11*XX1_p + p21*XX2_p) - sig*(I1 - (p11*PI1_p + p21*PI2_p))+ sig*RN;
f(end+1) = -PI1 + kapp*XX1 +bet*(p11*PI1_p + p21*PI2_p) +U; 
f(end+1) = -I1 + psi1*PI1 + psi_x*XX1 + IB;

f(end+1) = -XX2 + (p12*XX1_p + p22*XX2_p) - sig*(I2 - (p12*PI1_p + p22*PI2_p))+ sig*RN;
f(end+1) = -PI2 + kapp*XX2 +bet*(p12*PI1_p + p22*PI2_p) +U; 
f(end+1) = -I2 + psi2*PI2 + psi_x*XX2 + IB;

f(end+1) = IB_p - rho_i*IB;
f(end+1) = RN_p - rho_r*RN;
f(end+1) = U_p - rho_u*U;

%Check Computation of Steady-State Numerically
fnum = double(subs(f, [Y X YP XP], [ss, ss]));
disp('Checking steady-state equations:')
disp(fnum);

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
