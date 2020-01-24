% materials 14
% Goals:
% 1.) maybe a last attempt at getting rid of the overshooting
% 2.) a correct way of dealing with endogenous states
% 3.) a correct way of doing the projection facility
% 4.) a correction for PLM that leads to new fa & fb
% 23 Jan 2020
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


%% Do IRFs for extensions - same code as before
command_IRFs_many_learning
