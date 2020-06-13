% generate motivation plots
% 9 Feb 2020

% clearvars
% close all
clc

date_today = strrep(datestr(today, 'yyyy-mm-dd'), '-','_');

% % Add all the relevant paths and grab the codename
this_code = mfilename;
[current_dir, basepath, BC_researchpath,toolpath,export_figpath,figpath,tablepath,datapath] = add_paths;
[fs, lw] = plot_configs;

print_figs = 0;
do_infl_exp =0; % this one take a couple of seconds extra for some reason
close_em =0;

%% unemployment rate

% series_id = 'UNRATE';
% observation_start = '2010-01-01';
% observation_end   = datestr(today,'yyyy-mm-dd');
% [output] = getFredData(series_id, observation_start, observation_end);
% urate = output.Data(:,2);
% time_urate = output.Data(:,1);

% figure
% set(gcf,'color','w'); % sets white background color
% set(gcf, 'Position', get(0, 'Screensize')); % sets the figure fullscreen
% plot(time_urate,urate, 'linewidth',lw)
% ax = gca; % current axes
% ax.FontSize = fs;
% set(gca,'TickLabelInterpreter', 'latex');
% grid on
% grid minor
% datetick('x','yyyy-mm', 'keeplimits')
% if print_figs ==1
%     figname = ['urate_', date_today]
%     cd(figpath)
%     export_fig(figname)
%     cd(current_dir)
%     close
% end


% this thing doesn't seem to work with figures
file='create_publishable_code_w_figures.m';
options = struct('format','pdf','outputDir','/Users/lauragati/Dropbox/BC_Research/next/code/publish_code');
docpath = publish(file,options)