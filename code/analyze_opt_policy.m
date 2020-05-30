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

%%
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
pibvals = linspace(-0.05, 0.05, n); % most values in simulations are between these
k1mean  = mean(k1grid)*ones(1,n);
pmean = zeros(1,n);
smean = zeros(1,n);

% approximate policy function
ppi = csapi({k1grid,pgrid,sgrid,sgrid,sgrid,sgrid},it);

% hold every at their means, just move k1
i_k = fnval(ppi,{k1vals,pmean,smean,smean,smean,smean});
plot(k1vals,i_k(:,1,1,1,1,1))

i_k = fnval(ppi,[k1vals;pmean;smean;smean;smean;smean]);
plot(k1vals,i_k) % this gives you the same thing - sure b/c the interpretation is that you have a history of zeros, but the gain moves
% almost like a partial derivative, this is.

% do the same for pibar, but keep k1 at its mean
i_p = fnval(ppi,[k1mean;pibvals;smean;smean;smean;smean]);
plot(pibvals,i_p)

create_simple_plot_x(k1vals,i_k,'i(k^{-1})','',[this_code, '_ik'],print_figs)
create_simple_plot_x(pibvals,i_p,'i(\bar{\pi})','',[this_code, '_ip'],print_figs)

% Let's do a getting unanchored scenario towards inflation and deflation
pibvals = linspace(0,0.5,n);
i_unanch1 = fnval(ppi,[k1vals;pibvals;smean;smean;smean;smean]);
create_simple_plot_x(pibvals,i_unanch1,'i','',[this_code, '_unanch_inf'],print_figs)

pibvals = linspace(0,-0.5,n);
i_unanch2 = fnval(ppi,[k1vals;pibvals;smean;smean;smean;smean]);
create_simple_plot_x(pibvals,i_unanch2,'i','',[this_code, '_unanch_def'],print_figs)
