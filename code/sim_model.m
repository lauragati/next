% sim_model > simulate data from model solution
% very strongly inspired by Ryan's sim_dat.m
% 12 sept 2019


function [xsim, ysim] = sim_model(gx,hx,eta,nper,ndrop,shock,x0)
ny = size(gx,1);
nx = size(hx,1);

ysim = zeros(ny,nper);
xsim = zeros(nx,nper);

%Initial state
if nargin>6
    xsim(:,1) = x0;
else
    xsim(:,1) = zeros(nx,1);
end

for j = 1:nper-1
    ysim(:,j) = gx*xsim(:,j);
    xsim(:,j+1) = hx*xsim(:,j)+eta*shock(:,j+1);
end
ysim(:,j+1) = gx*xsim(:,j+1);

xsim = xsim(:,ndrop+1:end);
ysim = ysim(:,ndrop+1:end);
