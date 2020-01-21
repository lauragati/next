function [IR, iry, irx]=ir_RE_Markov_switching(gx,hx,r,x0,T)
%IR=ir(gx,hx,x0,T) computes T-period 
%impulse responses of the 
%vector [y x], whose law of 
%motion is:
% x(t+1) = hx x(t)
% y(t) = gx x2(t)
% with initial condition x0
% Inputs: gx, hx, x0, and 
%T (optional, default value 10)
% Adapted from 
%(c) Stephanie Schmitt-Grohe and Martin Uribe, August 18, 1993. 
% by Laura Gati 21 Jan 2020.

if nargin < 5
T =10;
end 
ny = size(gx,1)/2;
gx1 = gx(1:ny,:);
gx2 = gx(ny+1:end,:);

x0=x0(:);
pd=length(x0);
% MX=[gx;eye(pd)];
IR=[];
x=x0;
for t=1:T
    if r(t)==1 % active regime
        MX=[gx1;eye(pd)];
    elseif r(t)==2 % passive regime
        MX=[gx2;eye(pd)];
    end
IR(t,:)=(MX*x)';
x = hx * x;
end
% ny = size(gx,1);
iry = IR(:,1:ny);
irx = IR(:,ny+1:end); 