% MODEL_LINEAR - Linearized version of the SOE without a binding collateral
% constraint.  Model is equivalent to SGU 2003 with (a) a slightly
% different specification of adjustment costs (b) a different definition of
% investment
%
% usage:
%
% [fy, fx, fyp, fxp] = model_linear(param, xidx, yidx, xpidx, ypidx)
%
% where
%
% param = a parameters object, created in parmeters.m
%
% NOTES: This program, which generates the model matrices in cannonical form,
%        requires the MATLAB symbolic toolbox! The algorithms used to solve 
%        the model, however, are numerical, and require no additonal m-files
%        beyond those that appear in this directory.
%
% Code by Ryan Chahrour, Boston College, 2014
 

function [mod,param,set] = model(param,set)


%Name of text files
mod.fname   = 'model_prog.m';
mod.ss_call = 'model_ss.m';

%Steady State

[ss param] = model_ss(param,set);

%Declare parameters symbols: parameters are values to be estimated
param_list = fieldnames(param);
syms(param_list{1:end});
for j = 1:length(param_list)
    eval(['PARAM(j) = ',param_list{j} ,';']);
end
PARAM

%Declare setting symbols: symbols are values that are "fixed"
set_list = fieldnames(set);
syms(set_list{1:end});
for j = 1:length(set_list)
    eval(['SET(j) = ',set_list{j} ,';']);
end
SET


%Declare Needed Symbols
syms A A_p K K_p  GAM GAM_p
syms C C_p H H_p W W_p R R_p I I_p

%Declare X and Y vectors
X  = [K   GAM ];
XP = [K_p GAM_p ];

Y  = [C   H   W   R   I];
YP = [C_p H_p W_p R_p I_p] ;


%Loop creates indexes to track variables (for figures, etc)
ny = length(Y);
nx = length(X);
for j=1:ny; 
    vn = char(Y(j));
    eval([lower(vn), '_idx = j;'])
end
for j=1:nx; 
    vn = char(X(j));
    eval([lower(vn), '_idx = ny+j;'])
end

%Model Equations
f(1) = C-W/chi;
f(2) = 1-bet*C/C_p*GAM_p^(-1/(1-alph))*(R_p+1-del);
f(3) = R - alph*GAM*(K/H)^(alph-1);
f(4)  = W - (1-alph)*GAM^(alph/(alph-1))*(K/H)^alph;
f(5) = (1-del)*K + I*GAM^(1/(1-alph)) - K_p*GAM^(1/(1-alph));
f(6) = C+I- GAM^(alph/(alph-1))*K^alph*H^(1-alph);
f(7) = log(GAM_p/gam);

%Log-linear approx (Pure linear if log_var = [])
xlog = 1:length(X);
ylog = 1:length(Y); ylog(1:2) = []; %H in levels for policy
log_var = [X(xlog) Y(ylog) XP(xlog) YP(ylog)];



mod.f = subs(f, log_var, exp(log_var));
mod.X = X;
mod.XP = XP;
mod.Y = Y;
mod.YP = YP;
mod.PARAM = PARAM;
mod.param = param;
mod.SET = SET;
mod.set = set;
mod.adiff = 0; %Include anaylytical derivatives?
mod.xlog = xlog;
mod.ylog = ylog;


%Standard Errors
mod.shck = sym(zeros(length(X),4));


%Measurement Error (parameters are std-devs in param.m file)
mod.me = [];

%Derivatives using numerical toolbox
mod = anal_deriv(mod);

%Save index variables for use in main_prog
!rm -f v_idx.mat
save v_idx *_idx


%*************************************************************************
% MAKE_PRIME - For any symoblic vector, return a vector with symoblic
% elements of the same names with a _p suffix.
%
%
% usage (by example)
%
% out = make_prime([X1 X2])
%
% out =
%       [X1_p X2_p]
%**************************************************************************
function Vp = make_prime(V)

Vp = sym('Vp');

for j = 1:length(V)
    Vp(j) = sym([char(V(j)), '_p']);
end


