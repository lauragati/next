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
print_figs        = 0;
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

% return

% figure
% set(gcf,'color','w'); % sets white background color
% set(gcf, 'Position', get(0, 'Screensize')); % sets the figure fullscreen
% plot(time, infl, 'linewidth', lw); hold on
% plot(time, spf1(1:end-1), 'linewidth', lw);
% plot(time, spf10(1:end-1), 'linewidth', lw);
% ax = gca; % current axes
% ax.FontSize = fs;
% datetick('x','yyyy', 'keeplimits')
% % The next three lines force the figure to start where the data starts
% xaxislimits= get(gca,'XLim');
% xaxislimits(1) = time(1);
% set(gca, 'XLim', xaxislimits);
% set(gca,'TickLabelInterpreter', 'latex');
% grid on
% grid minor
% legend('CPI inflation (yoy, \%)', '1-year ahead forecast, SPF (annual, \%)', '10-year expected average, SPF (annual, \%)', 'location', 'southoutside', 'interpreter', 'latex')
% legend('boxoff')


% % return
% figure
% set(gcf,'color','w'); % sets white background color
% set(gcf, 'Position', get(0, 'Screensize')); % sets the figure fullscreen
% plot(time, fe, 'linewidth', lw); hold on
% % plot(x,zeros(1,T), 'k--', 'linewidth',lw)
% ax = gca; % current axes
% ax.FontSize = fs;
% datetick('x','yyyy', 'keeplimits')
% % The next three lines force the figure to start where the data starts
% xaxislimits= get(gca,'XLim');
% xaxislimits(1) = time(1);
% set(gca, 'XLim', xaxislimits);
% set(gca,'TickLabelInterpreter', 'latex');
% grid on
% grid minor
% legend('Forecast errors (yoy, \%)', 'location', 'southoutside', 'interpreter', 'latex')
% legend('boxoff')
% 
% 
% figure
% histogram(fe)

%% Regress LR-E on fe
% first half of sample
mdl = fitlm(fe(1:T/2),spf10(2:T/2+1))

% second half of sample
mdl = fitlm(fe(T/2+1:end),spf10(T/2+2:end))

% omit first fraction of sample
mdl = fitlm(fe(28+1:end),spf10(28+2:end))

% time(28+1) '1999-Q1'

% Let's implement a rolling regression
n_windows = 3;
length_window = T/n_windows;
idx_window = 1:T/n_windows:T;
% add the last
idx_window = [idx_window, T];

bet_win = zeros(n_windows,1);
for i=1:n_windows
    X = fe(idx_window(i)+1:idx_window(i+1));
    Y = spf10(idx_window(i)+2:idx_window(i+1)+1);
    bet_win(i) = (X'*X) \ (X'*Y);
end
% not providing great evidence here


% Try changes
% full sample
mdl = fitlm(fe,spf10(2:end)-spf10(1:end-1))

% first half of sample
mdl = fitlm(fe(1:T/2),spf10(2:T/2+1)-spf10(2-1:T/2+1-1) )

% second half of sample
mdl = fitlm(fe(T/2+1:end),spf10(T/2+2:end)-spf10(T/2+2-1:end-1))
% no they're not good

% return

%% NY FED SCE

xlsx_file = '/Users/lauragati/Dropbox/BC_Research/next/data/raw/NY_Fed_SCE/edited.xlsx';

sheet = 'expectations';
[num,txt,raw] = xlsread(xlsx_file, sheet);
% Median 3-year ahead inflation expectations
sce_3median = num(:,end);

sheet = 'uncertainty';
[num,txt,raw] = xlsread(xlsx_file, sheet);
% Median uncertainty around 3-year ahead
sce_unc_median = num(:,end);

sheet = 'distr';
% % of responses in the ranges (-Inf, 0)	[0,1)	[1,2)	[2,3)	[3,4)	[4, Inf)
[num,txt,raw] = xlsread(xlsx_file, sheet);
sce_distr = num(:,8:end);

dates_sce = num(:,1);
dates_sce = num2str(dates_sce);
dates_sce = [dates_sce(:,1:4), repmat('-',size(dates_sce,1),1), dates_sce(:,5:end)];
dates_sce = datenum(dates_sce);
dates_sce_str = datestr(dates_sce, 'yyyy-mm');

figure
set(gcf,'color','w'); % sets white background color
set(gcf, 'Position', get(0, 'Screensize')); % sets the figure fullscreen
plot(dates_sce, sce_3median, 'linewidth', lw); hold on
% plot(x,zeros(1,T), 'k--', 'linewidth',lw)
ax = gca; % current axes
ax.FontSize = fs;
datetick('x','yyyy', 'keeplimits')
% The next three lines force the figure to start where the data starts
xaxislimits= get(gca,'XLim');
xaxislimits(1) = dates_sce(1);
set(gca, 'XLim', xaxislimits);
set(gca,'TickLabelInterpreter', 'latex');
grid on
grid minor
legend('NY Fed SCE median 3-year ahead inflation expectations (yoy, \%)', 'location', 'southoutside', 'interpreter', 'latex')
legend('boxoff')

figure
set(gcf,'color','w'); % sets white background color
set(gcf, 'Position', get(0, 'Screensize')); % sets the figure fullscreen
plot(dates_sce, sce_unc_median, 'linewidth', lw); hold on
% plot(x,zeros(1,T), 'k--', 'linewidth',lw)
ax = gca; % current axes
ax.FontSize = fs;
datetick('x','yyyy', 'keeplimits')
% The next three lines force the figure to start where the data starts
xaxislimits= get(gca,'XLim');
xaxislimits(1) = dates_sce(1);
set(gca, 'XLim', xaxislimits);
set(gca,'TickLabelInterpreter', 'latex');
grid on
grid minor
legend('NY Fed SCE median 3-year ahead uncertainty', 'location', 'southoutside', 'interpreter', 'latex')
legend('boxoff')

figure
set(gcf,'color','w'); % sets white background color
set(gcf, 'Position', get(0, 'Screensize')); % sets the figure fullscreen
for d=1:size(sce_distr,2)
    h(d) = plot(dates_sce, sce_distr(:,d), 'linewidth', lw); hold on
end
ax = gca; % current axes
ax.FontSize = fs;
datetick('x','yyyy', 'keeplimits')
% The next three lines force the figure to start where the data starts
xaxislimits= get(gca,'XLim');
xaxislimits(1) = dates_sce(1);
set(gca, 'XLim', xaxislimits);
set(gca,'TickLabelInterpreter', 'latex');
grid on
grid minor
% title('NY Fed SCE 3-year ahead bins (yoy, \%)', 'interpreter', 'latex')
legend(h, '(-$\infty$,0)', '[0,1)','[1,2)', '[2,3)', '[3,4)', '(4, $\infty$)','location', 'southoutside', 'interpreter', 'latex', 'NumColumns',6)
legend('boxoff')


% Replot the distribution for the 0,1 and 3,4 bins
figure
set(gcf,'color','w'); % sets white background color
set(gcf, 'Position', get(0, 'Screensize')); % sets the figure fullscreen
hlow = plot(dates_sce, sce_distr(:,2), 'linewidth', lw); hold on
hhigh = plot(dates_sce, sce_distr(:,end-1), 'linewidth', lw);
ax = gca; % current axes
ax.FontSize = fs;
datetick('x','yyyy', 'keeplimits')
% The next three lines force the figure to start where the data starts
xaxislimits= get(gca,'XLim');
xaxislimits(1) = dates_sce(1);
set(gca, 'XLim', xaxislimits);
set(gca,'TickLabelInterpreter', 'latex');
grid on
grid minor
% title('NY Fed SCE 3-year ahead bins (yoy, \%)', 'interpreter', 'latex')
legend([hlow, hhigh], '[0,1)', '[3,4)', 'location', 'southoutside', 'interpreter', 'latex', 'NumColumns',2)
legend('boxoff')


%% Livingston


xlsx_file = '/Users/lauragati/Dropbox/BC_Research/next/data/raw/Livingston/medians.xlsx';

sheet = 'CPI';
[num,txt,raw] = xlsread(xlsx_file, sheet);
% Median 3-year ahead inflation expectations
liv_10y = num(:,end);
liv_10y = liv_10y(~isnan(liv_10y));


xlsx_file = '/Users/lauragati/Dropbox/BC_Research/next/data/raw/Livingston/Dispersion1.xlsx';

sheet = 'CPI10Y';
[num,txt,raw] = xlsread(xlsx_file, sheet);
% Median 3-year ahead inflation expectations
liv_10y_IQR = num(:,end);
liv_10y_IQR = liv_10y_IQR(~isnan(liv_10y_IQR));

start_liv = datenum('1991-Q2', 'yyyy-qq');
end_liv = datenum('2020-Q2', 'yyyy-qq');

time_liv = linspace(start_liv, end_liv,numel(liv_10y));
datestr(time_liv, 'yyyy-qq'); % it has small errors

figure
set(gcf,'color','w'); % sets white background color
set(gcf, 'Position', get(0, 'Screensize')); % sets the figure fullscreen
plot(time_liv, liv_10y, 'linewidth', lw); hold on
% plot(x,zeros(1,T), 'k--', 'linewidth',lw)
ax = gca; % current axes
ax.FontSize = fs;
datetick('x','yyyy', 'keeplimits')
% % The next three lines force the figure to start where the data starts
% xaxislimits= get(gca,'XLim');
% xaxislimits(1) = dates_sce(1);
% set(gca, 'XLim', xaxislimits);
set(gca,'TickLabelInterpreter', 'latex');
grid on
grid minor
legend('Livingston median 10-year ahead CPI inflation expectations (\%)', 'location', 'southoutside', 'interpreter', 'latex')
legend('boxoff')

figure
set(gcf,'color','w'); % sets white background color
set(gcf, 'Position', get(0, 'Screensize')); % sets the figure fullscreen
plot(time_liv, liv_10y_IQR, 'linewidth', lw); hold on
% plot(x,zeros(1,T), 'k--', 'linewidth',lw)
ax = gca; % current axes
ax.FontSize = fs;
datetick('x','yyyy', 'keeplimits')
% % The next three lines force the figure to start where the data starts
% xaxislimits= get(gca,'XLim');
% xaxislimits(1) = dates_sce(1);
% set(gca, 'XLim', xaxislimits);
set(gca,'TickLabelInterpreter', 'latex');
grid on
grid minor
legend('Livingston 10-year ahead CPI inflation expectations, interquartile range', 'location', 'southoutside', 'interpreter', 'latex')
legend('boxoff')

%% inflation expectations, market based
url = 'https://fred.stlouisfed.org/'; % frequency is daily!
c = fred(url);
series = 'T10YIE';
d = fetch(c,series);
time_breakeven10  = d.Data(:,1);
breakeven10 = d.Data(:,2);

series = 'T20YIEM';
d = fetch(c,series);
time_breakeven20  = d.Data(:,1);
breakeven20 = d.Data(:,2);

series = 'T30YIEM';
d = fetch(c,series);
time_breakeven30  = d.Data(:,1);
breakeven30 = d.Data(:,2);

close(c)



%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%% For prezi & paper
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

if print_figs==1
    close all
else
    disp('not displaying prezi plots')
    return
end

%% Forecast errors in the SPF data

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
xlabel('Forecast error time series, \%', 'interpreter', 'latex', 'fontsize', fs)

% legend('Forecast errors (yoy, \%)', 'location', 'southoutside', 'interpreter', 'latex')
% legend('boxoff')

figname = ['fe_SPF_', figspecs];

if print_figs ==1
    disp(figname)
    cd(figpath)
    export_fig(figname)
    cd(current_dir)
    close
end


figure
set(gcf,'color','w'); % sets white background color
set(gcf, 'Position', get(0, 'Screensize')); % sets the figure fullscreen
histogram(fe)
ax = gca; % current axes
ax.FontSize = fs;
set(gca,'TickLabelInterpreter', 'latex');
grid on
grid minor
xlabel('Forecast error, \%', 'interpreter', 'latex', 'fontsize', fs)


figname = ['fe_SPF_hist_', figspecs];

if print_figs ==1
    disp(figname)
    cd(figpath)
    export_fig(figname)
    cd(current_dir)
    close
end

% Compute mean and median fe 
mean(fe)
sort_fe = sort(fe);
if mod(T,2)==1
    median_fe = sort_fe(ceil(T/2));
else
    median_fe = (sort_fe(T/2)+sort_fe(T/2+1))/2
end

% Compute implied gain time series for the  "complete Materials 44 candidate" (11 Oct 2020)
fegrid = [-4,-3,0,3,4];
% map to ndim_simplex
x = cell(1,1);
x{1} = fegrid;
alph_draft_sept2125 = [0.8161,    0.6133,    0,    0.3342,    0.4452]'; % "complete Materials 44 candidate"
k_ts = ndim_simplex_eval(x,fe',alph_draft_sept2125);

% Compute mean and median gain 
mean(k_ts)
sort_k_ts = sort(k_ts);
if mod(T,2)==1
    median_k_ts = sort_k_ts(ceil(T/2));
else
    median_k_ts = (sort_k_ts(T/2)+sort_k_ts(T/2+1))/2
end

% Plot ts and histogram of gain
figure
set(gcf,'color','w'); % sets white background color
set(gcf, 'Position', get(0, 'Screensize')); % sets the figure fullscreen
plot(time, k_ts, 'linewidth', lw); hold on
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
xlabel('Gain time series', 'interpreter', 'latex', 'fontsize', fs)

% legend('Forecast errors (yoy, \%)', 'location', 'southoutside', 'interpreter', 'latex')
% legend('boxoff')

figname = ['gain_SPF_', figspecs];

if print_figs ==1
    disp(figname)
    cd(figpath)
    export_fig(figname)
    cd(current_dir)
    close
end


figure
set(gcf,'color','w'); % sets white background color
set(gcf, 'Position', get(0, 'Screensize')); % sets the figure fullscreen
histogram(k_ts)
ax = gca; % current axes
ax.FontSize = fs;
set(gca,'TickLabelInterpreter', 'latex');
grid on
grid minor
xlabel('Gain', 'interpreter', 'latex', 'fontsize', fs)


figname = ['gain_SPF_hist_', figspecs];

if print_figs ==1
    disp(figname)
    cd(figpath)
    export_fig(figname)
    cd(current_dir)
    close
end

% Compute implies time series of g(.)*fe
% Plot ts and histogram of gain
figure
set(gcf,'color','w'); % sets white background color
set(gcf, 'Position', get(0, 'Screensize')); % sets the figure fullscreen
plot(time, fe.*k_ts, 'linewidth', lw); hold on
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
xlabel('$\mathbf{g}(\cdot) f_{t|t-1}$ time series (annualized percentage points)', 'interpreter', 'latex', 'fontsize', fs)

% legend('Forecast errors (yoy, \%)', 'location', 'southoutside', 'interpreter', 'latex')
% legend('boxoff')

figname = ['gdot_fe_SPF_', figspecs];

if print_figs ==1
    disp(figname)
    cd(figpath)
    export_fig(figname)
    cd(current_dir)
    close
end

% do a single plot that synthetizes fe, gain and g(.)*fe
% Plot ts and histogram of gain
figure
set(gcf,'color','w'); % sets white background color
set(gcf, 'Position', get(0, 'Screensize')); % sets the figure fullscreen
yyaxis left
plot(time, fe, 'linewidth', lw, 'linestyle', ':'); hold on
plot(time, fe.*k_ts, 'linewidth', lw,'linestyle', '-');
ylabel('Annualized percentage points', 'interpreter', 'latex')
yyaxis right
plot(time, k_ts, 'linewidth', lw,'linestyle', '-.'); 
ax = gca; % current axes
ax.FontSize = fs*4/5;
datetick('x','yyyy', 'keeplimits')
% The next three lines force the figure to start where the data starts
xaxislimits= get(gca,'XLim');
xaxislimits(1) = time(1);
set(gca, 'XLim', xaxislimits);
set(gca,'TickLabelInterpreter', 'latex');
grid on
grid minor
% xlabel('Gain time series', 'interpreter', 'latex', 'fontsize', fs)

legend('Forecast errors', 'Change in long-run expectations', 'Gain', 'location', 'southoutside', 'interpreter', 'latex', 'NumColumns', 3, 'fontsize', fs*4/5)
legend('boxoff')

figname = ['fe_gain_gdot_fe_SPF_', figspecs];

if print_figs ==1
    disp(figname)
    cd(figpath)
    export_fig(figname)
    cd(current_dir)
    close
end

% Compute mean and median g(.)*fe 
gdotfe = fe.*k_ts;
sort_gdotfe = sort(gdotfe);
mean_gdotfe = mean(gdotfe) % -0.0659
if mod(T,2)==1
    median_gdotfe = sort_gdotfe(ceil(T/2));
else
    median_gdotfe = (sort_gdotfe(T/2)+sort_gdotfe(T/2+1))/2 % -0.0078
end

min(gdotfe) % -2.2433
max(gdotfe) % 0.7458



%% The regression plot using CPI inflation and SPF 1- and 10-year ahead

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


figname = ['regression_plot_', figspecs];

if print_figs ==1
    disp(figname)
    cd(figpath)
    export_fig(figname)
    cd(current_dir)
    close
end

%% Plot Livingston and SPF 10 y-ahead next to each other
figure
set(gcf,'color','w'); % sets white background color
set(gcf, 'Position', get(0, 'Screensize')); % sets the figure fullscreen
hliv = plot(time_liv, liv_10y, 'linewidth', lw); hold on
hspf = plot(time, spf10(1:end-1), 'linewidth', lw);
% plot(x,zeros(1,T), 'k--', 'linewidth',lw)
ax = gca; % current axes
ax.FontSize = fs;
datetick('x','yyyy', 'keeplimits')
% % The next three lines force the figure to start where the data starts
% xaxislimits= get(gca,'XLim');
% xaxislimits(1) = dates_sce(1);
% set(gca, 'XLim', xaxislimits);
set(gca,'TickLabelInterpreter', 'latex');
grid on
grid minor
% title('Median 10-year ahead CPI inflation expectations (annual, \%)', 'interpreter', 'latex')
legend([hliv, hspf],'Livingston', 'Survey of Professional Forecasters (SPF)', 'location', 'southoutside', 'interpreter', 'latex')
legend('boxoff')

figname = ['epi_in_data_', figspecs];

if print_figs ==1
    disp(figname)
    cd(figpath)
    export_fig(figname)
    cd(current_dir)
    close
end

%% Plot 3 breakeven inflation series next to each other
figure
set(gcf,'color','w'); % sets white background color
set(gcf, 'Position', get(0, 'Screensize')); % sets the figure fullscreen
hbe10 = plot(time_breakeven10, breakeven10, 'linewidth', lw);hold on
hbe20 = plot(time_breakeven20, breakeven20, 'linewidth', lw);
hbe30 = plot(time_breakeven30, breakeven30, 'linewidth', lw);
% plot(x,zeros(1,T), 'k--', 'linewidth',lw)
ax = gca; % current axes
ax.FontSize = fs;
datetick('x','yyyy', 'keeplimits')
% The next three lines force the figure to start where the data starts
xaxislimits= get(gca,'XLim');
xaxislimits(1) = dates_sce(1);
set(gca, 'XLim', xaxislimits);
set(gca,'TickLabelInterpreter', 'latex');
grid on
grid minor
% title('Median 10-year ahead CPI inflation expectations (annual, \%)', 'interpreter', 'latex')
legend([hbe10, hbe20,hbe30],'10-year', '20-year','30-year', 'location', 'southoutside', 'interpreter', 'latex','NumColumns',3)
legend('boxoff')

figname = ['epi_be_in_data_', figspecs];

if print_figs ==1
    disp(figname)
    cd(figpath)
    export_fig(figname)
    cd(current_dir)
    close
end

%% Plot NY FED SCE distribution (3-year ahead infl-E)
% Replot the distribution for the 0,1 and 3,4 bins
figure
set(gcf,'color','w'); % sets white background color
set(gcf, 'Position', get(0, 'Screensize')); % sets the figure fullscreen
hlow = plot(dates_sce, sce_distr(:,2), 'linewidth', lw); hold on
hhigh = plot(dates_sce, sce_distr(:,end-1), 'linewidth', lw);
ax = gca; % current axes
ax.FontSize = fs;
datetick('x','yyyy', 'keeplimits')
% The next three lines force the figure to start where the data starts
xaxislimits= get(gca,'XLim');
xaxislimits(1) = dates_sce(1);
set(gca, 'XLim', xaxislimits);
set(gca,'TickLabelInterpreter', 'latex');
grid on
grid minor
% title('NY Fed SCE 3-year ahead bins (yoy, \%)', 'interpreter', 'latex')
legend([hlow, hhigh], '[0,1)', '[3,4)', 'location', 'southoutside', 'interpreter', 'latex', 'NumColumns',2)
legend('boxoff')

figname = ['SCE_distrib_topbottom_', figspecs];
if print_figs ==1
    disp(figname)
    cd(figpath)
    export_fig(figname)
    cd(current_dir)
    close
end


%% Livingston IQR (10-year ahead)

figure
set(gcf,'color','w'); % sets white background color
set(gcf, 'Position', get(0, 'Screensize')); % sets the figure fullscreen
plot(time_liv, liv_10y_IQR, 'linewidth', lw); hold on
% plot(x,zeros(1,T), 'k--', 'linewidth',lw)
ax = gca; % current axes
ax.FontSize = fs;
datetick('x','yyyy', 'keeplimits')
% % The next three lines force the figure to start where the data starts
% xaxislimits= get(gca,'XLim');
% xaxislimits(1) = dates_sce(1);
% set(gca, 'XLim', xaxislimits);
set(gca,'TickLabelInterpreter', 'latex');
grid on
grid minor
% legend('Livingston 10-year ahead CPI inflation expectations, interquartile range', 'location', 'southoutside', 'interpreter', 'latex')
% legend('boxoff')

figname = ['Livingston_IQR_', figspecs];
if print_figs ==1
    disp(figname)
    cd(figpath)
    export_fig(figname)
    cd(current_dir)
    close
end