%**************************************************************
% MAIN_PROG - Solves the neoclassical model with a random walk 
% TFP solved in class in EC860. Compute impulstylee responses, and 
% plots them.
%
% Code by Ryan Chahrour, Boston College, 2012
%**************************************************************

nper =10;

disp('SOLVING LINEAR MODEL ...')
load model_object

[param,set] = parameters;
set.logvars = [modl.ylog,modl.xlog];
[~,param,set] = model_ss(parameters,set);

%Compute the first-order coefficiencients of the model
paramv = struct2array(param);
setv = struct2array(set);
[f, fx, fy, fxp, fyp, eta]=model_prog(paramv,setv);
ss = model_ss(param,set);

%Check steady-state
if max(abs(f))>1e-8
    disp('Warning: Non-zero steady-state');
    pause
end
[gx,hx]=gx_hx_alt(fy,fx,fyp,fxp);


%Finish HERE!