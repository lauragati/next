%% Coax the solver to get to the right answer: for 1D case only

 p = genpath(pwd)

nsearch = 100;
filename = 'acf_sim_univariate_data_21_Jun_2020'; % simulated data, nfe = 6. full Om
load([filename, '.mat'])
alph_true = acf_outputs{8};
nfe=6;
k1min = 0;
k1max= 1;
femax = 5;
femin = -femax;

% Uniform random starting values
rng('default')
b=1; a=0;
ALPH0 = a + (b-a).*rand(nfe,nsearch);
alph_opt = ones(nfe,nsearch);
resnorm = zeros(1,nsearch);
res = zeros(45,nsearch);
Om_opt = zeros(45,nsearch);
flag = zeros(1,nsearch);

tic
parfor i=1:nsearch
    alph0 = ALPH0(:,i);
    [alph_opt(:,i), resnorm(i),res(:,i), Om_opt(:,i),flag(i)] = fun_GMM_LOMgain_univariate(acf_outputs, nfe, k1min, k1max, femin, femax, alph0);
end
toc

[min_resnorm,min_idx]= min(resnorm);

[alph_true, alph_opt(:,min_idx)]