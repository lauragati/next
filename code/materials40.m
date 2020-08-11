% materials40
% still: overall identification in the estimation. Now I'm not sure it's a
% bug. It could be. The good news is that there seems to be a global min in
% the loss around the truth. I need to figure out why it isn't getting
% there.
% 5 August 2020

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


save_stuff=0;



%%

command_GMM_LOMgain_univariate

command_acf_sim_data_univariate

command_check_simulation_approx