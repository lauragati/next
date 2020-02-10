% materials16

% Goals:
% 1.) Prepare simulated loss plots for Clough prezi
% 2.) Replicate DW plots for Clough paper (Draft 1) - generate paper
% settings

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

%% Plot simulated loss for various internal learning parameters
plot_sim_loss

%% Do IRFs, simulation and gains for anchoring model, CEMP or CUSUM criterion
command_IRFs_anchoring

%% Do fmincon for psi_pi in the anchoring model, now also RE
grid_search

%% Replicate DW plots urate, FFR, inflation and Epi
create_motivation_plots