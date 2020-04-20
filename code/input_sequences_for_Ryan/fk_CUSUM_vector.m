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

% estimate of FEV
om = omt_1 + kap*kt_1^(-1)*(f*f' - omt_1);

%% Updating step - this is the issue

% my old approach: pre-29 Jan 2020
% thet = thett_1 + mean(mean(kap*kt_1^(-1)*(om^(-1)*(f*f')- thett_1))); % 

% for scalar CUSUM: try srqt of the square (actually works for vector case too)
% thet = thett_1 + mean(mean(kap*kt_1^(-1)*(om^(-1)*(f)- thett_1)));


% Ryan, 29 Jan 2020
% thet = thett_1 + kap*kt_1^(-1)*(det(om^(-1)*(f*f'))- thett_1); % % det is tiny!!

% based on Lütkepohl, "Intro to Multiple Time Series Analysis", p. 161
tauhat = f'*om^(-1)*f; %<--- fave

ny = size(f,1);
tautilde = f'*om^(-1)*f/ny;
thet = thett_1 + kap*kt_1^(-1)*(tauhat - thett_1); % this gives you a scalar yo







% an approach I wanna avoid (29 Jan 2020)
% % a projection facility..? --> om is exploding and need to deal with this!
% if max(abs(eig(om))) >1
%     om = omt_1;
%     thet = thett_1;
% end

%% Compare to threshold
I = thet <= thettilde;
k = I.*(kt_1+1)+(1-I).*gbar^(-1);
