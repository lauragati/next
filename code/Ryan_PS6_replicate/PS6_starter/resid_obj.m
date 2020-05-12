%RESID_OBJ - Compute the residual of five of the model equations for given
%policy functions Gpol,hpol,Vpol,Epol,Hpol.

function [out] = resid_obj(coeff,pw,Xgrid,XXgrid,epsg,param,set,MM)

%Unpacks parameters for you.
param_unpacker

%Number of shocks...
neps = length(pw);
ns   = length(XXgrid);

%Complete here
out = [];