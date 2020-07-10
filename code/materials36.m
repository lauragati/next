% materials36
% Continue to estimate the 1D approximated anchoring function, with
% autocovariogram and different seeds for simulated data and estimation
% 24 June 2020

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

do21 = 0;
do22 = 0;
do23 = 0;
do24 = 0;
do25 = 1;

save_stuff=0;

%%

command_GMM_LOMgain_univariate