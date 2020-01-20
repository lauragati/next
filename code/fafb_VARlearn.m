function [fa, fb] = fafb_VARlearn(param, a, b, z)
% a version of LH expectations for "VAR learning"ар la Slobodyan and
% Wouters
% 20 Jan 2020

bet = param.bet;  
alph = param.alph;
ny = size(b,2);

mu = (eye(ny)-b)^(-1)*a;
fb = 1/(1-bet)*mu + (eye(ny)-bet*b)^(-1)*b^2*(z-mu);
fa = 1/(1-alph*bet)*mu + (eye(ny)-alph*bet*b)^(-1)*b^2*(z-mu);