% Smooth anchoring function for the case that only pi is learned
% Still isn't working great.
% 19 March 2020
function [k, g_pi, g_pibar] = fk_smooth_pi_only(param,fe, kt_1)
% d=param.d;
% c=param.c;
rho_k = param.rho_k;
gam_k = param.gam_k;

% Alternative 1:
% k1 = kt_1^(-1) + c + d*(fe^2);
% k = k1^(-1);

% Alternative 2:
% k = kt_1 + d*1/(fe^2);

% % Alternative 3: % used to be preferred, but decreases monotonically <---------
% d = 10; c=0;
% k = kt_1 + 1/((d*fe)^2);
% g_1_pi = -2*(d*fe)^(-2)*fe^(-1);
% g_1_pibar = 2*(d*fe)^(-2)*fe^(-1);
% g_pi =2*d*fe;
% g_pibar = -2*d*fe;

% % Alternative 4 % This ain't so bad in the cross-section, but
% % individual gain series can go negative, so that's not cool.
% d=1; c=0.5;
% k = kt_1 + d - c*(fe^2);

% % <----
% % Alternative 5: % <-- This one works nicely! Try to go with this.
% d = 10; c=0.01;
% k = kt_1 + 1/((d*fe)^2) - c*kt_1;

% Alternative 6 = compatible with target criterion (k_t = g(fe))
% Empirically, doesn't work at all because the gain is specified in levels.
% k = 1/((d*fe)^2);

% <----
% Alternative 7 % This one looks like constant gain learning, so that's
% also considerable. It never goes negative, it fluctuates up and down
% around some mean. So should follow up on this one too.
% This is the only one that works with fsolve.
% rho=0.9; %0.5 or 0.9
% gam=0.001; %0.01 or 0.001
k1 = rho_k * kt_1.^(-1) + gam_k*(fe).^2;
k=1./k1;
g_pi = 2*gam_k*fe;
g_pibar = -2*gam_k*fe;

% % Alternative 8 - to make Alt 7 decrease in general
% rho=0.5; %0.5
% gam=0.01; %0.01
% k1 = (rho * kt_1^(-1) + gam*(fe)^2)/t;
% k=1/k1;



%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% % try the following decreasing gain scheme
% k = kt_1 +1;
% 
% % try the following constant gain scheme
% k = kt_1;