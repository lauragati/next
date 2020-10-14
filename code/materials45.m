% materials45
% Correcting Sept 21/25 draft
% 8 Oct 2020

clearvars
close all
clc

% Add all the relevant paths and grab the codename
this_code = mfilename;
[current_dir, basepath, BC_researchpath,toolpath,export_figpath,figpath,tablepath,datapath, inputsRyan_path] = add_paths;
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


save_stuff=0;

return


%% Estimation and fake data generation

command_GMM_LOMgain_univariate

command_acf_sim_data_univariate

coax_solver_server % to run the same estimation from different starting points

%% Calibration of sigmas and alphas

command_sigmas

%% Anchoring in the data

command_anchoring_in_data

command_anchoring_in_data_individual

%% Optimal policy

command_pea_approx_main

command_vfi_approx

%% Analyze optimal policy

compare_value_pea_results_approx

analyze_opt_policy

%% Plot CB loss, find optimal Taylor-rule coefficient and compute welfare

plot_sim_loss_approx

grid_search_approx

command_welfare

%% IRFs

command_IRFs_approx_pretty

