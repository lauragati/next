%% Attempt to estimate CEMP parameters with Metropolis-Hastings
clearvars
% COPY OF LAST SECTION OF MATERIALS3.M, INTENDED FOR BEING RUN ON SERVER
disp('--------- Bayesian estimation of CEMP for server ---------')

load('cemp_data.mat')
% The data has this structure:
% Y = [pi, spf1,spf2,liv1,liv2]';

%non-nan sample:
Y = Y(:,182:end);


T = size(Y,2);
l = size(Y,1);
k_nonlin = 2; % number of nonlinear states
k_lin = 3; % number of linear states
ne = 2; % number of state innovations
k = k_nonlin+k_lin;
m = 1000; % number of particles (with 2500 particles, it takes 20 sec, with 1000, 7 sec.)

ub = [8, 1,    0.2,  0.9,  0.99, 0.99, 0.5, 0.5, 0.5, 0.5, 0.5, 0.5, 0.5  ];
lb = [1, 0.001,0.05, 0.05, 0.01, 0.01, 0.01,0.01,0.01,0.01,0.01,0.01,0.01 ];

[param0, pmean, pstd, pshape] = parameters_cemp;

%Generate an MCMC chain
step = (ub-lb)/200;
% The bigger the stepsize, the more rejections you have.

%inital point: using outcome of MLE for initial params
disp('--------- Initial evaluation of likelihood ---------')
tic
p0 = marginalized_particlefilter_cemp(Y,param0,m);
toc
lnprior0 =0;

return
nonflat_prior=1;
if nonflat_prior ==1
    % If you choose a nonflat prior, specify prior shape, prior mean and prior
    % std.dev. All of these need to be (nparam x 1).
    % 1= beta; 2=gamma; 3 = normal; 4 = inv wishart; 5 = uniform.
    lnprior0 = evaluate_prior_materials3(param0,pshape,pmean,pstd);
end
%Start chain
tic
% chain length:
chaint = 1000;%100000; One chain element takes about 20 sec (with 2000 particles). So 100,000 should take 24 days!! wtf! So 1000 length should take 6 hours.
% 1000 chaint and 1000 particles should be 2 hours. (Actually took 3 hrs.)
% number of acceptances:
accept = 0;
% chain of params:
pchain = zeros(chaint,length(param0));
priorchain = zeros(chaint,1);
pchain(1,:) = param0;
% chain of posterior probability as a function of chain of params
p = zeros(1,chaint);
p(1) = p0 + lnprior0;
priorchain(1) = lnprior0;
% (the last two are both initialized at output MLE)
for j = 1:chaint
    [pchain(j+1,:),p(j+1), accept_paramp, priorchain(j+1)] = mcmc_cemp(pchain(j,:),Y,p(j),step,lb,ub,pshape,pmean,pstd,m);
    accept = accept+accept_paramp;
    if mod(j,1000)==0
        disp(['iteration ', num2str(j/1000), ' thousand.']);
    end
end
accept_rate = accept/chaint
toc

% Bayesian estimate
meanp = mean(pchain);

savefilename = 'cemp_bayesian_estimation_outputs.mat';
save(savefilename, 'p','pchain')
disp('--------- Done. ---------')
