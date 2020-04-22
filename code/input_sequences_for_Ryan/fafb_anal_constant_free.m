function [fa, fb] = fafb_anal_constant_free(param, a, b, s, hx)
% with 'free' I mean that it's independent of n
% 19 Oct 2019
% Note: This code simply computes LH-expectations blindly - you need to
% specify elsewhere (in the model A-matrices) whether this
% learning-expectation is used for each of the observables. 

bet = param.bet;  
alph = param.alph;

nx = size(hx,1);

fa = a/(1-alph*bet) + b*(eye(nx)-alph*bet*hx)^(-1)*s;
fb = a/(1-bet) + b*(eye(nx)-bet*hx)^(-1)*s;