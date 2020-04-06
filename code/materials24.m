% materials24
% 2 April 2020

% Goals:
% 1.) Figure out a way to guess a sequence of sthg and make it compatible
% with the model via optimization (target crit will be a special case)

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

skip_old_stuff = 1;
%% 1) Simulate given a sequence - optimize over that sequence to satisfy model

command_sim_given_seq
% this needs to be corrected, there is an fsolve way to do it conceptually
% better

%% 2) Value function iteration to solve for optimal i-sequence
command_valfun_iter