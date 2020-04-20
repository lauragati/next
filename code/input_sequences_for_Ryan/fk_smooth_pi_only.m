% Smooth anchoring function for the case that only pi is learned
% Still isn't working great.
% 19 March 2020
function [k, g_pi, g_pibar] = fk_smooth_pi_only(param,fe, kt_1)
d=param.d;
c=param.c;


% Alternative 1:
% k1 = kt_1^(-1) + c + d*(fe^2);
% k = k1^(-1);
% % k=k1;
% Alternative 2:
% k = kt_1 + d*1/(fe^2);

% Alternative 3: % best so far <---------
k = kt_1 + 1/((d*fe)^2);
g_1_pi = -2*(d*fe)^(-2)*fe^(-1);
g_1_pibar = 2*(d*fe)^(-2)*fe^(-1);
g_pi =2*d*fe;
g_pibar = -2*d*fe;

% % % Alternative 4
% k = kt_1 + d - c*(fe^2);
% g_pi = -2*c*fe;
% % Alternative 5: 
% k = kt_1 + 1/((d*fe)^2) - c;
% g_pi = -2*(d*fe)^(-2)*fe^(-1);
% Alternative 6 = compatible with target criterion (k_t = g(fe))
% k = 1/((d*fe)^2);



% % try the following decreasing gain scheme
% k = kt_1 +1;
% 
% % try the following constant gain scheme
% k = kt_1;
