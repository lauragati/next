% MODEL_SS - Return the steady state of the model computed analytically
%
% usage:
% 
% [ss, parameters] =model_ss(param)

%START_ADD
function [XX, PI, W, H, Q, RP, D, C, GAMY, GAMH, GAMC, GAMQ, GAMD, GAMYP, YC, YCP, GAMHP, HC, HCP,GAM, S, HL, CL, QL, DL, GAMYL, YCL, YCLL, GAMHL, HCL, HCLL, EPS, RPL,XX_p, PI_p, W_p, H_p, Q_p, RP_p, D_p, C_p, GAMY_p, GAMH_p, GAMC_p, GAMQ_p, GAMD_p, GAMYP_p, YC_p, YCP_p, GAMHP_p, HC_p, HCP_p,GAM_p, S_p, HL_p, CL_p, QL_p, DL_p, GAMYL_p, YCL_p, YCLL_p, GAMHL_p, HCL_p, HCLL_p, EPS_p, RPL_p param set] = model_ss(param,set)
%Assign parameter values to named variables.
rho = param(1);
phi = param(2);
rhoa = param(3);
alph = param(4);
sigz = param(5);
sige = param(6);
%Assign set values to named variables.

bet = set(1);
gam = set(2);
eta = set(3);
pis = set(4);
mu = set(5);
hbar = set(6);
rho_r = set(7);
rbar = set(8);
lam = set(9);
%END_ADD

%Use closed form expressions for the ss values.
w = (eta-1)/eta;
h = (w/(1-rho/gam))^(mu/(mu+1));
c = h;
x = c*(1-rho/gam);
R = pis*gam/bet;
d = c/eta;
q = d*bet/(1-bet);
s = rho*c;

%Update the param structure
set(8) = R;
set(5) = mu;


%This part you have to do by hand :-(
XX=x;
PI=pis;
W=w;
H=h;
Q=q;
RP=R;
D=d; 
C=c;
GAMY=gam;
GAMH=1;
GAMC=gam; 
GAMQ=gam;
GAMD=gam;
GAMYP=gam;
YC=1;
YCP=1;
GAMHP=1;
HC=1;
HCP=1;
GAM=gam;
S=s;
HL=h;
CL=c;
QL=q;
DL=d;
GAMYL=gam;
YCL=1;
YCLL=1;
GAMHL=1;
HCL=1;
HCLL=1;
EPS=1;
RPL=R;

XX_p=x;
PI_p=pis;
W_p=w;
H_p=h;
Q_p=q;
RP_p=R;
D_p=d; 
C_p=c;
GAMY_p=gam;
GAMH_p=1;
GAMC_p=gam; 
GAMQ_p=gam;
GAMD_p=gam;
GAMYP_p=gam;
YC_p=1;
YCP_p=1;
GAMHP_p=1;
HC_p=1;
HCP_p=1;
GAM_p=gam;
S_p=s;
HL_p=h;
CL_p=c;
QL_p=q;
DL_p=d;
GAMYL_p=gam;
YCL_p=1;
YCLL_p=1;
GAMHL_p=1;
HCL_p=1;
HCLL_p=1;
EPS_p=1;
RPL_p=R;
-1.000000e+00
-1.000000e+00
-1.000000e+00
-1.000000e+00
-1.000000e+00
-1.000000e+00
-1.000000e+00
-1.000000e+00
-1.000000e+00
-1.000000e+00
-1.000000e+00
-1.000000e+00
-1.000000e+00