% Builds on main_prog in the starter kit to solve a 2-dimensional
% stochastic value function iteration problem using interpolation (for the endogenous state) and
% discretization for the exogenous one.

clearvars
close all
clc

% Grab the codename
this_code = mfilename;
% Add all the relevant paths
current_dir = pwd;
PS6_starter_path = [current_dir, '/PS6_starter'];
cd ../.. % go up 2 levels
basepath = pwd;
cd .. % go up another level to BC_Research
BC_researchpath = pwd;
toolpath = [BC_researchpath '/matlab_toolbox'];
export_figpath = [toolpath '/Export_Fig'];
figpath = [basepath '/figures'];
tablepath = [basepath '/tables'];
datapath = [basepath '/data'];
tryouts_path = [toolpath '/tryouts'];
maincode_path = [basepath '/code'];


cd(current_dir)

addpath(basepath)
addpath(toolpath)
addpath(export_figpath)
addpath(figpath)
addpath(datapath)
addpath(tryouts_path)
addpath(PS6_starter_path)
addpath(maincode_path)

todays_date = strrep(datestr(today), '-','_');

print_figs=0;
skip = 1;

%**************************************************************************
[param,set] = parameters;
param_unpack

ndim = 2;
nk = 25;
ng = 5;

%**************************************************************************
%SOLVE LINEAR MODEL
%**************************************************************************
solve_linear

%**************************************************************************
% 1a) Simulate the linearized economy
T=5000;
ndrop=0;
rng(0)
eta = 1; % I think there's only the tech shock
[yxsim, shock] = sim_dat(gx,hx,eta,T,ndrop);
ksim = yxsim(k_idx,:);
sd_k = sqrt(var(ksim));

%**************************************************************************
% 1b) Create a 25-point grid for capital and 5-point grid for tech
% growthrate.
ks = ss(k_idx);
kmax = log(ks) + 3*sd_k;
kmin = log(ks) - 3*sd_k;
nk = 25; % no of capital gridpoints
kgrid = linspace(kmin,kmax,nk);
ng=5; % number of tech gripoints
[~,ggrid,pg] = GH_Quadrature(ng,1,1); % for some reason need to treat variance as 1
ggrid = ggrid*sige; % now multiply grid with variance
ggrid = ggrid(end:-1:1)' + log(gam); % flip the grid and add log stst of tech grt
pg = pg(end:-1:1); % flip the weights too
disp('Confirming std of shocks:')
disp(['siga: ' num2str(sqrt((ggrid.^2)*pg))]); % I think the aha-effect here is that sig(a) on the grid is close to sige in the linear model.

%**************************************************************************
% 1c) Using ndgrid, generate a 2x125 gridmatrix (X). Use the policy
% functions from the linearized econ to generate labor policy values Hi for
% each point on the grid. Confirm that the average h-policy value is close
% to the ss-value of labor.
[X1,X2] = ndgrid(ggrid,kgrid);
[row, col] = size(X1);
X = zeros(2,row*col);
idx=0;
for i=1:row
    for j=1:col
        idx=idx+1;
        X(:,idx) = [X2(i,j);X1(i,j)]; % put kgrid on top, ggrid below
    end
end
% X equals Ryan's XXgrid in the sols
% Still need to demean to compare with hours
gs=ss(gam_idx);
XXdemean = X - repmat_col([ks;gs],length(X));

[~,m] = size(X); % # gridpoints (= Ryan's 'ns')

% Initial policy using policy functions of linearized econ (gx summarizes those)
Y = zeros(ny,m);
for i=1:m
    Y(:,i) = gx*X(:,i);
    %      Y(:,i) = gx*XXdemean(:,i); % this ain't cool either
end
H0 = Y(h_idx,:);
hs = ss(h_idx);
hs - mean(H0) % that doesn't look good. Maybe it isn't in levels? ignore it for now.
% well...

%**************************************************************************
% 1d) Using initial policy values, generate a set of coeffs a_i such that
% Htilde(a_i,X) approximates the data.

pp = csapi({ggrid,kgrid},reshape(H0,size(X1)));
a =pp.coefs;

%**************************************************************************
% 2a) Residual function. Call it once

res0 = res_obj(pp,X,kgrid,ggrid,pg,param,set);
% holy shit it works - although the values are crap

% handle
obj = @(coeffs)res_obj(a,X,kgrid,ggrid,pg,param,set);

%**************************************************************************
% 3) Solve

options = optimset('fsolve');
options = optimset(options, 'display', 'iter');

% POL_final = fsolve(obj,pp,options); % aha, shit, fsolve doesn't handle
% the darn pp object



function res = res_obj(coeffs,Xgrid,kgrid,ggrid,pg,param,set)
param_unpack
ng=length(ggrid);
nk=length(kgrid);
% 1. H(a,X) at each gridpoint X_i
H = fnval(coeffs,{ggrid,kgrid});
H = reshape(H,1,ng*nk);
% 2. Use eqs. 1 and 2 to generate kp
K = exp(Xgrid(1,:));
GAM = exp(Xgrid(2,:));
C = (1-alph)/chi.*GAM.^(alph/(alph-1)).*(K./H).^alph;
KP = GAM.^(alph/(alph-1)).*K.^alph.*H.^(1-alph) - C +(1-del).*K.*GAM.^(alph/(alph-1));
% 3. For each KP(X_i), generate 5 next period grid points XP_ij according
% to grid values of GAM.
[X1,X2] = ndgrid(ggrid,KP);
[row, col] = size(X1);
XP = zeros(2,row*col);
idx=0;
for i=1:row
    for j=1:col
        idx=idx+1;
        XP(:,idx) = [X2(i,j);X1(i,j)]; % put kgrid on top, ggrid below
    end
end
% 4. Evaluate policy function at the points XP. Use period t+1 versions of
% eqs 1 and 2 to compute the terms on the RHS of eq 3. (rhs_i).
HP = fnval(coeffs,{ggrid,log(KP)});
% [h1,h2] = size(HP);
% HP =reshape(HP,1,h1*h2);
CP = (1-alph)/chi.*GAM.^(alph/(alph-1)).*(KP./HP).^alph;
rhs=bet./CP.*GAM.^(1/(1-alph)).* ((KP./HP).^(alph-1) +1-del );
Erhs = pg'*rhs;
% 5. Compute residuals
res = 1./C - Erhs;
end

