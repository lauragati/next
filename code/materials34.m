% materials34
% Split up elements of materials33 to continue to estimate the approximated anchoring function
% 17 June 2020

clearvars
close all
clc

% Add all the relevant paths and grab the codename
this_code = mfilename;
[current_dir, basepath, BC_researchpath,toolpath,export_figpath,figpath,tablepath,datapath, inputsRyan_path] = add_paths;
todays_date = strrep(datestr(today), '-','_');
nowstr = strrep(strrep(strrep(datestr(now), '-','_'), ' ', '_'), ':', '_');

% Variable stuff ---
print_figs        = 1;
stop_before_plots = 0;
skip_old_plots    = 0;
output_table = print_figs;

skip = 1;
[fs, lw] = plot_configs;
redo_data_load_and_bootstrap = 0;
datestr(now)


%% Get the data, filter them, generate data moments and bootstrap to get the weighting matrix
if skip==0
command_acf_data
end
%% Estimate given a dataset

command_GMM_LOMgain

%%  Estimate univariate anchoring function
comman_GMM_LOMgain_univariate