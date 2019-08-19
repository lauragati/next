function [Rhat, Rhat_alt] = gelman_rubin(X)
% X is a set of simulated samples for one specific parameter x
% Rhat is the square root of Rhat in Gelman-Rubin 1992, p. 461.

[n,m] = size(X);
% n is the length of the simulated sample (initially simulate 2n and discard n)
% m is the number of simulated samples (for CEMP I think m=5)

meanx_i = mean(X,1); % within-sample means (m,1)
meanx = mean(mean(X)); % the overall mean of x (1,1)
B = n/(m-1)* sum((meanx_i - meanx).^2); % (1,1)

s_i = var(X,1); % within-sample variances
W = mean(s_i);

muhat = meanx;

sighat = (n-1)/n * W + 1/n * B;

% Vhat = sighat + B/(m*n); 
Vhat = sighat + B/(m)*n; 

vars_i = var(s_i);
VC1 = cov(s_i,meanx_i.^2);
cov1 = VC1(1,2);
VC2 = cov(s_i,meanx_i);
cov2 = VC2(1,2);
varVhat = ((n-1)/n)^2 * 1/m * vars_i + ((m+1)/(m*n))^2 * 2/(m-1)*B^2 + 2*(m+1)*(n-1)/(m*n^2) * n/m* (cov1 - 2*muhat*cov2);
df = 2*Vhat^2 / varVhat;

Rhat = sqrt(Vhat/W*df/(df-2));

Rhat_alt = sqrt(sighat/W); % based on https://astrostatistics.psu.edu/RLectures/diagnosticsMCMC.pdf (17 July 2019) 
% stata blog agress with this alternative thing
% also Patrick Lam's cool lecture notes on convergence in MCMC (saved as patrick_lam_convergence.pdf in literature)
