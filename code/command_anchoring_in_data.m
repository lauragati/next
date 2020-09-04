% command_anchoring_in_data.m
% Look at evidence of unanchoring in data (following Stephen Terry)
% 28 August 2020

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
datestr(now)


save_stuff=0;

%% Data (in logs)

% SPF Survey of professional forecasters: CPI inflation expectation 1- and
% 10year ahead (annual rates, %, quarterly frequency)
xlsx_file = '/Users/lauragati/Dropbox/BC_Research/next/data/raw/SPF/Inflation.xlsx';
[num,txt,raw] = xlsread(xlsx_file);
spf1  = num(:,4);
spf10 = num(:,5);

begin_i = find(isnan(spf10),1, 'last')+1;
spf1  = num(begin_i:end,4);
spf10 = num(begin_i:end,5);

time_spf = string(num(begin_i:end,[1,2]));
time_spf = strcat(time_spf(:,1),'q',time_spf(:,2));
time_spf = datenum(time_spf,'yyyyqq');
timestr_spf = datestr(time_spf,'yyyy-qq');

% Read in CPI
% observation_start = '1991-10-01'; % LR-E is the binding constraint at the beginning.
observation_start = '1991-01-01'; % LR-E is the binding constraint at the beginning.
observation_end   = '2020-04-01'; % real GDP is the binding constraint at the end
units = 'lin';
frequency = 'q';
aggregation_method = 'avg';

% CPI inflation
[output1] = getFredData('CPIAUCSL', observation_start, observation_end, units, frequency, aggregation_method);
cpi = output1.Data(:,2);
%
infl_yoy = (cpi(5:end) - cpi(1:end-4))./cpi(1:end-4)*100;
infl_qoq = (cpi(2:end) - cpi(1:end-1))./cpi(1:end-1)*100; % -> first inflation obs is for 1992-Q1
% annualized q-o-q percent change (See Annualizing Data from Dallas Fed)
% infl = ((infl_qoq/100+1).^4 -1)*100;
% infl = infl_qoq;
infl = infl_yoy;

% Construct forecast errors
fe = infl - spf1(1:end-1);



% Thus the sample is in the end from 1992-Q1 to 2020-Q2
y = [infl, spf1(1:end-1), fe];
[T,n] = size(y);
time = time_spf(2:end);
timestr = datestr(time,'yyyy-qq');



%% Do some nice interesting plots
figspecs = [this_code, '_', nowstr];

figure
set(gcf,'color','w'); % sets white background color
set(gcf, 'Position', get(0, 'Screensize')); % sets the figure fullscreen
plot(time, spf10(1:end-1), 'linewidth', lw); 
ax = gca; % current axes
ax.FontSize = fs;
datetick('x','yyyy', 'keeplimits')
% The next three lines force the figure to start where the data starts
xaxislimits= get(gca,'XLim');
xaxislimits(1) = time(1);
set(gca, 'XLim', xaxislimits);
set(gca,'TickLabelInterpreter', 'latex');
grid on
grid minor
% ylab = ylabel('\%', 'interpreter', 'latex');
% % Rotate ylabel and put it on top
% ylab.Position
% ylab.Rotation = 0
% ylab.Position(1) = ylab.Position(1) - 500; % move left
% ylab.Position(2) = ylab.Position(2)*1.3; % move up

legend('10-year expected average, SPF (annual, \%)', 'location', 'southoutside', 'interpreter', 'latex')
legend('boxoff')

figname = ['lr_epi_in_data_', figspecs];

if print_figs ==1
    disp(figname)
    cd(figpath)
    export_fig(figname)
    cd(current_dir)
    close
end

% return

figure
set(gcf,'color','w'); % sets white background color
set(gcf, 'Position', get(0, 'Screensize')); % sets the figure fullscreen
plot(time, infl, 'linewidth', lw); hold on
plot(time, spf1(1:end-1), 'linewidth', lw); 
plot(time, spf10(1:end-1), 'linewidth', lw); 
ax = gca; % current axes
ax.FontSize = fs;
datetick('x','yyyy', 'keeplimits')
% The next three lines force the figure to start where the data starts
xaxislimits= get(gca,'XLim');
xaxislimits(1) = time(1);
set(gca, 'XLim', xaxislimits);
set(gca,'TickLabelInterpreter', 'latex');
grid on
grid minor
legend('CPI inflation (yoy, \%)', '1-year ahead forecast, SPF (annual, \%)', '10-year expected average, SPF (annual, \%)', 'location', 'southoutside', 'interpreter', 'latex')
legend('boxoff')

figname = ['epi_in_data_', figspecs];

if print_figs ==1
    disp(figname)
    cd(figpath)
    export_fig(figname)
    cd(current_dir)
    close
end

% return
figure
set(gcf,'color','w'); % sets white background color
set(gcf, 'Position', get(0, 'Screensize')); % sets the figure fullscreen
plot(time, fe, 'linewidth', lw); hold on
% plot(x,zeros(1,T), 'k--', 'linewidth',lw)
ax = gca; % current axes
ax.FontSize = fs;
datetick('x','yyyy', 'keeplimits')
% The next three lines force the figure to start where the data starts
xaxislimits= get(gca,'XLim');
xaxislimits(1) = time(1);
set(gca, 'XLim', xaxislimits);
set(gca,'TickLabelInterpreter', 'latex');
grid on
grid minor
legend('Forecast errors (yoy, \%)', 'location', 'southoutside', 'interpreter', 'latex')
legend('boxoff')

figname = ['fe_in_data_', figspecs];

if print_figs ==1
    disp(figname)
    cd(figpath)
    export_fig(figname)
    cd(current_dir)
    close
end

figure
histogram(fe)

%% Regress LR-E on fe
% first half of sample
mdl = fitlm(fe(1:T/2),spf10(2:T/2+1))

% second half of sample
mdl = fitlm(fe(T/2+1:end),spf10(T/2+2:end))





