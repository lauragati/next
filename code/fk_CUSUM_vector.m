function [k, om, thet] = fk_CUSUM_vector(param,kt_1,omt_1, thett_1, f)
% reworked version of fk_cusum
% accomodates vector learning
% output: 
% k = k
% om is the firms' estimate of the FEV of y (ny x ny)
% thet is the firms' estimate of the variance-scaled FEV of y (ny x ny) - I
% take the mean of this to be the scalar statistic, so thet becomes a
% scalar.
% 26 Jan 2020

gbar = param.gbar;
kap = param.kap; % 0.80; allows the test-statistic to be revised at a different rate than the estimate of the mean. 0 < kap < 1.
thettilde = param.thettilde; %1.60; the new thetbar. I just set both to match CEMP's criterion approx.

om = omt_1 + kap*kt_1^(-1)*(f*f' - omt_1);
thet = thett_1 + mean(mean(kap*kt_1^(-1)*(om^(-1)*(f*f')- thett_1)));

% % a projection facility..? --> om is exploding and need to deal with this!
% if max(abs(eig(om))) >1
%     om = omt_1;
%     thet = thett_1;
% end

I = thet <= thettilde;
k = I.*(kt_1+1)+(1-I).*gbar^(-1);
