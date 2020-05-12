function [f fx fy fxp fyp eta R set dgam_dtheta deta_dtheta dR_dtheta xlag ylag] = model_prog(param, set)

%Assign parameter values to named variables.
sige = param(1);

%Assign set values to named variables.
bet = set(1);
chi = set(2);
del = set(3);
alph = set(4);
phi = set(5);
gam = set(6);
kbar = set(7);

%BEGIN_EXTRACT_HERE

%Use closed form expressions for the ss values.
r = (1/(bet*gam^(1/(1-alph))))-1+del;
kh = (r/(alph*gam))^(1/(alph-1));
w = (1-alph)*gam^(alph/(alph-1))*(kh)^alph;
ik = (gam^(1/(1-alph))-(1-del))/(gam^(1/(1-alph)));
c = w/chi;
ck = gam^(alph/(alph-1))*(kh)^(alph-1)-ik;

k = ck^-1*c;
i = ik*k;
h = kh^-1*k;

%Put the ss values in a vector consistent with Y and X vectors in model.m
Yss  = [c h w r i];
Xss  = [k gam];
ss  = [Yss Xss];


%END_EXTRACT_HERE
%Compute Steady State
K= Xss(1);
GAM= Xss(2);
C= Yss(1);
H= Yss(2);
W= Yss(3);
R= Yss(4);
I= Yss(5);
K_p= Xss(1);
GAM_p= Xss(2);
C_p= Yss(1);
H_p= Yss(2);
W_p= Yss(3);
R_p= Yss(4);
I_p= Yss(5);

%Evaluate F.
f = [[C - W/chi, 1 - (C*bet*(R_p - del + 1))/(C_p*GAM_p^(1/(alph - 1))), R - GAM*alph*(K/H)^(alph - 1), W + GAM^(alph/(alph - 1))*(K/H)^alph*(alph - 1), I/GAM^(1/(alph - 1)) - K*(del - 1) - K_p/GAM^(1/(alph - 1)), C + I - GAM^(alph/(alph - 1))*H^(1 - alph)*K^alph, log(GAM_p/gam)]];
%Evaluate derivative expressions.
fx = [[0, 0]; [0, 0]; [-(GAM*K*alph*(K/H)^(alph - 2)*(alph - 1))/H, -GAM*alph*(K/H)^(alph - 1)]; [(GAM^(alph/(alph - 1))*K*alph*(K/H)^(alph - 1)*(alph - 1))/H, GAM*GAM^(alph/(alph - 1) - 1)*alph*(K/H)^alph]; [-K*(del - 1), (GAM*K_p)/(GAM^(1/(alph - 1) + 1)*(alph - 1)) - (GAM*I)/(GAM^(1/(alph - 1) + 1)*(alph - 1))]; [-GAM^(alph/(alph - 1))*H^(1 - alph)*K*K^(alph - 1)*alph, -(GAM*GAM^(alph/(alph - 1) - 1)*H^(1 - alph)*K^alph*alph)/(alph - 1)]; [0, 0]];
fy = [[1, 0, -W/chi, 0, 0]; [-(bet*(R_p - del + 1))/(C_p*GAM_p^(1/(alph - 1))), 0, 0, 0, 0]; [0, (GAM*K*alph*(K/H)^(alph - 2)*(alph - 1))/H^2, 0, R, 0]; [0, -(GAM^(alph/(alph - 1))*K*alph*(K/H)^(alph - 1)*(alph - 1))/H^2, W, 0, 0]; [0, 0, 0, 0, I/GAM^(1/(alph - 1))]; [1, (GAM^(alph/(alph - 1))*K^alph*(alph - 1))/H^alph, 0, 0, I]; [0, 0, 0, 0, 0]];
fxp = [[0, 0]; [0, (C*GAM_p*bet*(R_p - del + 1))/(C_p*GAM_p^(1/(alph - 1) + 1)*(alph - 1))]; [0, 0]; [0, 0]; [-K_p/GAM^(1/(alph - 1)), 0]; [0, 0]; [0, 1]];
fyp = [[0, 0, 0, 0, 0]; [(C*bet*(R_p - del + 1))/(C_p^2*GAM_p^(1/(alph - 1))), 0, 0, -(C*R_p*bet)/(C_p*GAM_p^(1/(alph - 1))), 0]; [0, 0, 0, 0, 0]; [0, 0, 0, 0, 0]; [0, 0, 0, 0, 0]; [0, 0, 0, 0, 0]; [0, 0, 0, 0, 0]];

eta = [[0, 0, 0, 0]; [0, 0, 0, 0]];
R = [];
