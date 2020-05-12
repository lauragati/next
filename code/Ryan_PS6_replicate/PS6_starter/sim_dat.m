%SIM_DAT - Simulate date from the model solution
%
%usage
%
%[yxsim, shock] = sim_dat(gx,hx,eta,nper,ndrop)

function [yxsim, shock] = sim_dat(gx,hx,eta,nper,ndrop,x0,dist,shock)
ny = size(gx,1);
nx = size(hx,1);
neps = size(eta,2);

if nargin<7 || strcmp(dist, 'normal') && nargin<8  %If nargin == 8, shocks are input arg
    shock = randn(neps,nper);
elseif strcmp(dist, 'tdist');
    df = 3;
    shock = trnd(df,neps,nper)/sqrt(df/(df-2));
elseif strcmp(dist, 'nct')
    shock = nct(neps,nper);
elseif strcmp(dist, 'sknrnd')
    shock = sknrnd(-8,neps,nper);
end


ysim = zeros(ny,nper);
xsim = zeros(nx,nper);

%Initial state
if nargin>5
    xsim(:,1) = x0;
end

for j = 1:nper-1
    ysim(:,j) = gx*xsim(:,j);
    xsim(:,j+1) = hx*xsim(:,j)+eta*shock(:,j+1);
end
ysim(:,j+1) = gx*xsim(:,j+1);

xsim = xsim(:,ndrop+1:end);
ysim = ysim(:,ndrop+1:end);
shock = shock(:,ndrop+1:end);

yxsim = [ysim;xsim];