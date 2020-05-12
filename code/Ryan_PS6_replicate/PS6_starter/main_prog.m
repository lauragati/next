%**************************************************************
% MAIN_PROG - Solves the neoclassical model with shocks to 
% TFP solved in class in ECON8860. This code will get you started!
%
% Code by Ryan Chahrour, Boston College, 2016
%**************************************************************

[param,set] = parameters;
param_unpack

ndim = 2;
nk = 25;
ng = 5;

%**************************************************************************
%SOLVE LINEAR MODEL
%**************************************************************************
solve_linear



%**************************************************************************
%QUADRATURE POINTS (use GH_Quadrature)
%**************************************************************************


%**************************************************************************
%GRID POINTS (use ndgrid)
%**************************************************************************

%THIS IS THE APPROXIMATION GRID
Xgrid{1} = [];
Xgrid{2} = [];


%THESE ARE THE EVALUATION POINTS (2-by-125 matrix)


%**************************************************************************
%INITIAL POLICY FUNCTIONS (use ndim_simplex)
%**************************************************************************

%Replace this line with code initializing hours at each grid point
H_init = ones(125,1);

%This generate the coefficients a_i
[POL_init,MM] = ndim_simplex(Xgrid,XXgrid,H_init);

%This evalues the policy function at random points (for example)
[H_vals] = ndim_simplex_eval(Xgrid,randn(2,100),POL_init)

disp(' ')
disp('Compare init policy and steady state')
[mean(POL_init), hbar]

%**************************************************************************
%SOLVE
%**************************************************************************
options = optimset('fsolve');
options = optimset(options, 'display', 'iter', 'tolfun', 1e-9, 'maxfunevals', 50000, 'tolx', 1e-10, 'jacobian', 'on', 'algorithm', 'trust-region-dogleg', 'derivativecheck ', 'off');

obj = @(coeff) resid_obj(coeff,pw,Xgrid,XXgrid,epsg,paramv,setv,MM);

disp(' ');
disp('Testing the objective using initial values.')
max(abs(obj(POL_init)))

options = optimset('fsolve');
options = optimset(options, 'display', 'iter');

POL_final = fsolve(obj,POL_init,options);


%**************************************************************************
%PLOT POLICY FUNCTION
%**************************************************************************



%**************************************************************************
%SIMULATE ECONOMY
%**************************************************************************


%My code for making nice simulation tables
disp('NONLINEAR SIMULATED MOMENTS:')
mom_tab_sim(longsim_nl,vvidx,vlist,prefs);


disp('NONLINEAR/LINEAR CORRELATIONS:')
diag(corr(longsim', longsim_nl'))