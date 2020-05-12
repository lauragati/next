%REPMAT_ROW - Reproduce a row vector n-times.  Faster than matlab's
%built in function repmat.m
%
% usage:
%
% out = repmat_row(x,n)
%
% where
% 
% x = row vector
% n = number of replications
%
% NO ERROR CHECKING!
function out = rempat_row(x,n)

out = x(ones(1,n),:);


