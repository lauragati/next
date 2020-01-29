% materials15

% Goals:
% 1.) understand difference in  behavior CEMP vs. CUSUM criterion
% 29 Jan 2020
clearvars
close all
clc

% Add all the relevant paths and grab the codename
this_code = mfilename;
[current_dir, basepath, BC_researchpath,toolpath,export_figpath,figpath,tablepath,datapath] = add_paths;

% Variable stuff ---
print_figs        = 0;
stop_before_plots = 0;
skip_old_plots    = 0;
output_table = print_figs;

skip_old_stuff = 1;

%% Do IRFs, simulation and gains for anchoring model, CEMP or CUSUM criterion
command_IRFs_anchoring

%% Do fmincon for psi_pi in the anchoring model
grid_search