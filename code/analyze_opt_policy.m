% analyze_opt_policy.m
% Load value function iteration and some parameterized expectations results
% to analyze the policy function
% 30 May 2020

clearvars
close all
clc

% Add all the relevant paths and grab the codename
this_code = mfilename;
[current_dir, basepath, BC_researchpath,toolpath,export_figpath,figpath,tablepath,datapath, inputsRyan_path] = add_paths;
todays_date = strrep(datestr(today), '-','_');

% Variable stuff ---
print_figs        = 0;
stop_before_plots = 0;
skip_old_plots    = 0;
output_table = print_figs;

skip = 1;

%% Comparative statics of policy wrt k1 and pibar
[param, set, param_names, param_values_str, param_titles] = parameters_next;
sig_r = param.sig_r;
sig_u = param.sig_u;

% Grids
ns = 2;
sgrid = linspace(-sig_r,sig_r,ns);

% value_output_name = 'value_outputs.mat';
% value_output_name = 'value_outputs_server.mat';
% value_output_name = 'value_outputs_server32_accelerated';
value_output_name = 'value_outputs_30_May_2020_10_42_12';
% value_output_name = 'value_outputs_approx13_Jun_2020_17_49_01';
load([value_output_name, '.mat'])

pp     = value_sols{1};
v      = value_sols{2};
it     = value_sols{3};
pibp   = value_sols{4};
k1p    = value_sols{5};
pgrid  = value_sols{6};
k1grid = value_sols{7};

n = 10;
k1vals  = linspace(min(k1grid),max(k1grid),n);
% pibvals = linspace(min(pgrid),max(pgrid),n);

% get some means for pibvals
% take the history of states from parametric expectations
% pea_output_name = 'pea_outputs_30_May_2020_10_18_28';
pea_output_name = 'pea_outputs_29_May_2020_14_42_56';
% load('inputs.mat')
load([pea_output_name, '.mat'])
e = output{1};
ysim7 = output{2};
k7    = output{3};
phi7  = output{4};
seq_opt = output{5};
i_pe = seq_opt(3,:);

% compile state vector
T=length(e)-2; % drop the first and last obs
k1sim = 1./k7(2:end-1);
pibsim = squeeze(phi7(1,1,2:end-1))';
max(pibsim) % 0.0156, 0.2136
mean(pibsim) % -0.0012, 0.1945
min(pibsim) % -0.0438, -0.0931
ssim = e(:,2:end-1);

pibvals = linspace(-0.05, 0.05, n); % most values in simulations are between these
k1mean  = mean(k1grid)*ones(1,n);
pmean = zeros(1,n);
smean = zeros(1,n);

% approximate policy function
ppi = csapi({k1grid,pgrid,sgrid,sgrid,sgrid,sgrid},it);

% % hold every at their means, just move k1
% i_k = fnval(ppi,{k1vals,pmean,smean,smean,smean,smean});
% plot(k1vals,i_k(:,1,1,1,1,1))

i_k = fnval(ppi,[k1vals;pmean;smean;smean;smean;smean]);
% plot(k1vals,i_k) % this gives you the same thing - sure b/c the interpretation is that you have a history of zeros, but the gain moves
% almost like a partial derivative, this is.

% do the same for pibar, keep k1 at its mean
i_p = fnval(ppi,[k1mean;pibvals;smean;smean;smean;smean]);
% plot(pibvals,i_p)
% legend('$ i(k^{-1}) $', 'interpreter', 'latex') % math mode needs latex interpreter
% title('Title', 'interpreter', 'latex') % std latex font needs latex
% interpreter
% title('\fontname{Courier}{Title}', 'interpreter', 'tex') % changing fonts
% requires tex interpreter


% create_simple_plot_x(k1vals,i_k,'$ i(k^{-1}) $','',[this_code, '_ik'],print_figs)
% create_simple_plot_x(pibvals,i_p,'$ i(\bar{\pi}) $','',[this_code, '_ip'],print_figs)
create_pretty_plot_x(k1vals,i_k,[this_code, '_ik', todays_date],print_figs)
create_pretty_plot_x(pibvals,i_p,[this_code, '_ip', todays_date],print_figs)

% % Let's do a getting unanchored scenario towards inflation and deflation
% pibvals = linspace(0,0.05,n);
% i_unanch1 = fnval(ppi,[k1vals;pibvals;smean;smean;smean;smean]);
% create_simple_plot_x(pibvals,i_unanch1,'$i$','',[this_code, '_unanch_inf'],print_figs)
%
% pibvals = linspace(0,-0.05,n);
% i_unanch2 = fnval(ppi,[k1vals;pibvals;smean;smean;smean;smean]);
% create_simple_plot_x(pibvals,i_unanch2,'$i$','',[this_code, '_unanch_def'],print_figs)

return
%% Regress i on pi and x
X = [ysim7(1,2:end-1); ysim7(2,2:end-1)]';
Y = i_pe';
bet_ols = (X'*X)^(-1)*(X'*Y);
yhat = X*bet_ols;
ybar = mean(Y);
% sum of squared residuals:
SSR = sum((Y-yhat).^2);
% total sum of squares:
SST = sum((Y-ybar).^2);
R2 = 1 - SSR/SST;

%% Check how k1 and pibar comove
% remove the initial, negative pibar value
X = [k1sim(2:end)', ssim(1,2:end)', ssim(3,2:end)'];
Y = abs(pibsim(2:end)');
bet_ols = (X'*X)^(-1)*(X'*Y);
yhat = X*bet_ols;
ybar = mean(Y);
% sum of squared residuals:
SSR = sum((Y-yhat).^2);
% total sum of squares:
SST = sum((Y-ybar).^2);
R2 = 1 - SSR/SST;

% explained sum of squares
SSE = sum((yhat-ybar).^2);
R2_2 = SSE/SST;

cov(k1sim',pibsim')
corr(k1sim',pibsim')

% if beta = 5, then for every unit change in k1, pibar changes by 5.
% so what's the average increase in the sample in k1?
dk1 = diff(k1sim);
dk1plus = dk1(dk1>0);
dk1minus = dk1(dk1<0);
mean_increase_k1 = mean(dk1plus);
mean_response_pibar = 5*mean_increase_k1;

pibvals = linspace(-mean_response_pibar,mean_response_pibar,n);
k1vals =  linspace(-mean_increase_k1,mean_increase_k1,n);
i_mean_unanchor = fnval(ppi,[k1vals;pibvals;smean;smean;smean;smean]);
