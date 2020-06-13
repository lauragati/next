% generate motivation plots
% 9 Feb 2020

clearvars
close all
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

series_id = 'UNRATE';
observation_start = '2010-01-01';
observation_end   = datestr(today,'yyyy-mm-dd');
[output] = getFredData(series_id, observation_start, observation_end);
urate = output.Data(:,2);
time_urate = output.Data(:,1);

figure
set(gcf,'color','w'); % sets white background color
set(gcf, 'Position', get(0, 'Screensize')); % sets the figure fullscreen
plot(time_urate,urate, 'linewidth',lw)
ax = gca; % current axes
ax.FontSize = fs;
set(gca,'TickLabelInterpreter', 'latex');
grid on
grid minor
datetick('x','yyyy-mm', 'keeplimits')
if print_figs ==1
    figname = ['urate_', date_today]
    cd(figpath)
    export_fig(figname)
    cd(current_dir)
    close
end

%% PCE inflation

series_id = 'PCEPI';
observation_start = '2009-01-01';
observation_end   = datestr(today,'yyyy-mm-dd');
[output] = getFredData(series_id, observation_start, observation_end);
pce = output.Data(:,2);
% %-change from a year ago
pi_pce = (pce(13:end) - pce(1:end-12))./pce(1:end-12)*100;
time_pce = output.Data(:,1);
time_pce =time_pce(13:end);

figure
set(gcf,'color','w'); % sets white background color
set(gcf, 'Position', get(0, 'Screensize')); % sets the figure fullscreen
plot(time_pce,pi_pce, 'linewidth',lw); hold on
plot(time_pce, 2*ones(length(pi_pce)), 'k--', 'linewidth',lw)
ax = gca; % current axes
ax.FontSize = fs;
set(gca,'TickLabelInterpreter', 'latex');
grid on
grid minor
datetick('x','yyyy-mm', 'keeplimits')
if print_figs ==1
    figname = ['pce_', date_today]
    cd(figpath)
    export_fig(figname)
    cd(current_dir)
    close
end

%% Fed funds rate target upper limit
series_id = 'DFEDTARU';
observation_start = '2010-01-01';
observation_end   = datestr(today-60,'yyyy-mm-dd'); % the -60 is a trick to make the x-axis not extend into 2022
[output] = getFredData(series_id, observation_start, observation_end);
ffr_t = output.Data(:,2);
time_fedfunds = output.Data(:,1);

figure
set(gcf,'color','w'); % sets white background color
set(gcf, 'Position', get(0, 'Screensize')); % sets the figure fullscreen
plot(time_fedfunds,ffr_t, 'linewidth',lw)
ax = gca; % current axes
ax.FontSize = fs;
set(gca,'TickLabelInterpreter', 'latex');
grid on
grid minor
datetick('x','yyyy-mm', 'keeplimits')
if print_figs ==1
    figname = ['frr_', date_today]
    cd(figpath)
    export_fig(figname)
    cd(current_dir)
    close
end

if do_infl_exp==1
    %% inflation expectations, market based
    url = 'https://fred.stlouisfed.org/';
    c = fred(url);
    series = 'T10YIE';
    d = fetch(c,series);
    time_epi  = d.Data(:,1);
    epi = d.Data(:,2);
    close(c)
    
    obs_start= datenum(observation_start, 'yyyy-mm-dd');
    % observation_end   = datestr(today,'yyyy-mm-dd');
    observation_end   = datestr(today-63,'yyyy-mm-dd');
    obs_end= datenum(observation_end, 'yyyy-mm-dd');
    
    ind_start = find(time_epi==obs_start);
    ind_end = find(time_epi==obs_end);
    
    epi = epi(ind_start:ind_end);
    time_epi = time_epi(ind_start:ind_end);
    
    
    figure
    set(gcf,'color','w'); % sets white background color
    set(gcf, 'Position', get(0, 'Screensize')); % sets the figure fullscreen
    plot(time_epi,epi, 'linewidth',lw); hold on
    plot(time_epi, 2*ones(length(epi)), 'k--', 'linewidth',lw)
    ax = gca; % current axes
    ax.FontSize = fs;
    set(gca,'TickLabelInterpreter', 'latex');
    grid on
    grid minor
    datetick('x','yyyy-mm', 'keeplimits')
    if print_figs ==1
        figname = ['epi10_', date_today]
        cd(figpath)
        export_fig(figname)
        cd(current_dir)
        close
    end
end
%%
if close_em==1
    close all
end
