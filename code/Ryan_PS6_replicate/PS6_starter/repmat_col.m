%REPMAT_COL - Reproduce a column vector n-times.  Faster than matlab's
%built in function repmat.m
%
% usage:
%
% out = repmat_col(x,n)
%
% where
% 
% x = column vector
% n = number of replications
%
% NO ERROR CHECKING!

function out = rempat_col(x,n)

out = x(:,ones(n,1));


