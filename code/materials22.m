% materials22
% 19 March 2020

% Goals:
% 1.) Do initial GMM estimation of the anchoring function in midsimple
% model - think I can take some shortcuts from general codes
% 2.) Do the numerical implementation of the target criterion in midsimple
% model - same shortcuts apply?

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