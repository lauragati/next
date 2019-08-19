% I'm taking the equations from my work on Lorenzoni 2009

function [fyn, fxn, fypn, fxpn] = model_preston(param)

%Steady State
[ss, param] = model_preston_ss(param);

%Declare parameters
% f       = param.f;
c       = param.C; % IB is on top, RN is below
alph    = param.alph;
bet     = param.bet;
sig     = param.sig;
kapp    = param.kapp;
psi_x   = param.psi_x;
psi_pi  = param.psi_pi;
w_small = param.w_small;

%Declare Needed Symbols
syms RN RN_p IB IB_p
syms PI PI_p XX XX_p I I_p

%Declare X and Y vectors
X  = [RN IB]; % I'm splitting the shock vector S into r^n and ibar.
XP = [RN_p IB_p];

Y  = [PI XX I  ];
YP = [PI_p XX_p I_p ] ;


%Model Equations (Preston, prior to eq. 13 & 14)
f(1) = -XX + XX_p - sig*(I - PI_p)+ RN;
f(2) = -PI + kapp*XX +bet*PI_p; 
f(3) = -I + psi_pi*PI + psi_x*XX + IB;
f(4) = IB_p - c(1,1)*IB;
f(5) = RN_p - c(2,2)*RN;

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
