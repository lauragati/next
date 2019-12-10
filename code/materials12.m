% materials 12
% Goals:
% Find the changes in the model that make it fit data
% 1.) Changes to the policy rule
% 2.) Changes to expectation formation
% 5 Dec 2019
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



%% Do IRFs for vector learning and EE learning
%  E(pi) in TR instead of pi also implemented in the same code
command_IRFs_many_learning

%% lagged inflation in Taylor rule

