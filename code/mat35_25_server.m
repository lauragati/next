% Do point 2.5) in materials35 faster on the server
% 27 June 2020

clearvars
close all
clc

% Add all the relevant paths and grab the codename
this_code = mfilename;
% Server-paths
current_dir = pwd;
cd ../ % go up 1 levels
basepath = pwd;
cd([basepath, '/materials35/'])
addpath([basepath, '/matlab_toolbox'])
addpath([basepath, '/next_toolbox'])


todays_date = strrep(datestr(today), '-','_');
nowstr = strrep(strrep(strrep(datestr(now), '-','_'), ' ', '_'), ':', '_');

% Variable stuff ---
print_figs        = 0;
stop_before_plots = 0;
skip_old_plots    = 0;
output_table = print_figs;

skip = 1;
[fs, lw] = plot_configs;
datestr(now)

nsearch=10;
disp(['Expected to take ', num2str(nsearch*30/60), ' minutes.'])
filename ='acf_data_11_Jun_2020'; % real data
% filename = 'acf_sim_univariate_data_21_Jun_2020'; % simulated data, nfe = 6. full Om
load([filename, '.mat'])
nfe=5;
k1min = 0;
k1max= 1;
femax = 3.5;
femin = -femax;

% Uniform random starting values
rng('default')
% b=1; a=0;
% ALPH0 = a + (b-a).*rand(nfe,nsearch);
ALPH0 = rand(nfe,nsearch);
alph_opt = ones(nfe,nsearch);
resnorm = zeros(1,nsearch);
res = zeros(45,nsearch);
Om_opt = zeros(45,nsearch);
flag = zeros(1,nsearch);

% just used to check against command_GMM_LOMgain_univariate.m
% ALPH0 =     [0.0674
%     0.0168
%     0
%     0.0168
%     0.0674]; % default*5


tic
for i=1:nsearch
    alph0 = ALPH0(:,i);
    [alph_opt(:,i), resnorm(i),res(:,i), Om_opt(:,i),flag(i)] = fun_GMM_LOMgain_univariate(acf_outputs, nfe, k1min, k1max, femin, femax, alph0);
end
toc

alph_converged = alph_opt(:,flag>0);
resnorm_converged = resnorm(flag>0);
[min_resnorm,min_idx]= min(resnorm_converged);
alph_best = alph_converged(:,min_idx);
disp('The best candidate and mean of candidates:')
[alph_best, mean(alph_converged,2)]

filename = [this_code, '_best_candidates_', todays_date];
outputs = {alph_opt, resnorm, Om_opt, flag, alph_converged, resnorm_converged};
save([filename,'.mat'],'loss_outputs')
disp(['Saving as ' filename])