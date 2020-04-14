% Smooth anchoring function for vector or scalar learning (the case that
% only pi is learned)
% 19 March 2020
function [k, g_pi, g_pibar] = fk_smooth(param,fe, kt_1)
d=param.d;
c=param.c;


% % this will be the k^(-1) formulation, doesn't work yet
% g = @(fe) d*fe'*fe - c;
% k1 = 1/kt_1 + g(fe);
% k = 1/k1;

% Alternative 3: % best so far <---------
k = kt_1 + 1/((d*fe')*(d*fe)); % this at least should implement the vector version


% these are already for the k^(-1) formulation:
g_pi =2*d*fe;
g_pibar = -2*d*fe;
