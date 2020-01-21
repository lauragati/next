% MODEL_SS - Return the steady state of the model (computed analytically)
%
% usage:
% 
% [ss, param] = model_ss(param)


function [ss,param] = model_NK_Markov_switchingTR_ss(param)

%Use closed form expressions for the ss values.
pi1 = 0; % shouldn't zero st.st. inflation imply pi = 1?
i1  = 0;
x1  = 0;

pi2 = 0; 
i2  = 0;
x2  = 0;

ib = 0;
rn = 0;
u  = 0;
%Put the ss values in a vector consistent with Y and X vectors in model.m
yy  = [pi1 x1 i1 pi2 x2 i2];
xx  = [rn ib u];
ss  = [yy xx];