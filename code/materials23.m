% materials23
% 26 March 2020

% Goals:
% 1.) Continure initial GMM estimation of the anchoring function in midsimple
% model with improvements
% 2.) Continue the numerical implementation of the target criterion in midsimple
% model, maybe try a value function iteration approach too.

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
%% 1) Estimation of anchoring function

command_GMM_anchoring_function

%% 2) Implementation of optimal policy

command_implement_target_criterion
