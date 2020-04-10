% materials25
% Prepare macro lunch of 15 April 2020
% --> get some kind of implementation of Taylor rule to work
% 1. Return to fsolve with the cleaned up learning-simulation code

clearvars
close all
clc

% Add all the relevant paths and grab the codename
this_code = mfilename;
[current_dir, basepath, BC_researchpath,toolpath,export_figpath,figpath,tablepath,datapath] = add_paths;
todays_date = strrep(datestr(today), '-','_');

% Variable stuff ---
print_figs        = 0;
stop_before_plots = 0;
skip_old_plots    = 0;
output_table = print_figs;

skip = 1;
%% parked stuff
if skip==0
    command_valfun_iter
    play_value_function_iter
end

%% 1. Return to fsolve
command_sim_given_seq_fsolve