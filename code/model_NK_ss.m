% MODEL_SS - Return the steady state of the model (computed analytically)
%
% usage:
% 
% [ss, param] = model_ss(param)


function [ss,param] = model_NK_ss(param)

%Use closed form expressions for the ss values.
pi = 0; % shouldn't zero st.st. inflation imply pi = 1?
i  = 0;
x  = 0;
ib = 0;
rn = 0;
u  = 0;
%Put the ss values in a vector consistent with Y and X vectors in model.m
yy  = [pi x i];
xx  = [rn ib u];
ss  = [yy xx];
