% materials46.m
% A code that adds standard errors to the estimation of alphas, Fig. 2 in
% June and October draft.
% 6 June 2021

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% For a good overview of the project codes, see materials45.m.
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
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


%% Create bootstrapped samples 
% actually, this has already been done for 10000 samples in
% command_acf_data.m, so just use the first 100.

%% Run estimation routine on each bootstrap sample
% Do on server

command_GMM_boot_server
% saved estimates are called alpha_hat_boot_date.mat (and should be 5x100)


%% Take sddev of estimated alpha_hat
% do locally

se = command_se

%% Generate the plot
% do locally
% rely on command_sigmas.m


