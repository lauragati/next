% materials25
% Prepare macro lunch of 15 April 2020
% --> get some kind of implementation of Taylor rule to work
% 1. Return to fsolve with the cleaned up learning-simulation code
% that doesn't work

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

[fs, lw, grey, silver, maroon, light_coral, light_salmon,dark_green, green, light_green, light_sky_blue, ...
    teal, purple, saddle_brown, light_brown, fs_pres,lw_pres] = plot_configs;

%%
% UMich Inflation Expectations
series_id = 'MICH';
observation_start = '1978-01-01';
observation_end   = '1990-01-01';
[output] = getFredData(series_id, observation_start, observation_end);
umich = output.Data(:,2);
time_umich = output.Data(:,1);


figure
set(gcf,'color','w'); % sets white background color
set(gcf, 'Position', get(0, 'Screensize')); % sets the figure fullscreen
plot(time_umich,umich, 'k', 'linewidth',lw_pres)
ax = gca; % current axes
ax.FontSize = fs_pres;
grid on
grid minor
datetick('x','yyyy')
if print_figs ==1
    figname = ['macrolunch_umich']
    cd(figpath)
    export_fig(figname)
    cd(current_dir)
    close
end
