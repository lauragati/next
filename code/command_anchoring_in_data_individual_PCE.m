% command_anchoring_in_data_individual_PCE.m
% Look at evidence of unanchoring in data from individual forecasters (following Jenny Tang)
% for PCE and PCE core
% 20 Oct 2020

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

%% Data

%% Read in PCE and PCE core
% observation_start = '2007-01-01'; % PCE SR-E and LR-E is the binding constraint at the beginning.
observation_start = '2006-04-01'; % need to take a little earlier to construct inflation and forecast errors
observation_end   = '2020-07-01'; 
units = 'lin';
frequency = 'q';
aggregation_method = 'avg';

% PCE inflation
[output1] = getFredData('PCEPI', observation_start, observation_end, units, frequency, aggregation_method);
pce = output1.Data(:,2);
%
infl_yoy = (pce(5:end) - pce(1:end-4))./pce(1:end-4)*100;
infl_qoq = (pce(2:end) - pce(1:end-1))./pce(1:end-1)*100; % -> first inflation obs is for 1992-Q1
% annualized q-o-q percent change (See Annualizing Data from Dallas Fed)
% infl = ((infl_qoq/100+1).^4 -1)*100;
% infl = infl_qoq;
infl_pce = infl_yoy;

% PCE core inflation
[output1] = getFredData('DPCCRV1Q225SBEA', '2007-01-01', observation_end, units, frequency, aggregation_method);
infl_pce_core = output1.Data(:,2);


%% Read in individual LR forecasts from SPF
% (annual rates, %, quarterly frequency)
% Ten-Year Annual-Average Headline PCE Inflation (Horizon: Year of Survey and Following Nine Years)
xlsx_file = '/Users/lauragati/Dropbox/BC_Research/next/data/raw/SPF/individual_PCE10.xlsx';
[num,txt,raw] = xlsread(xlsx_file);
start_i = 6104-1; % minus to get rid of top row in excel
year = num(start_i:end,1);
quarter = num(start_i:end,2);

time_spf = string([year,quarter]);
time_spf = strcat(time_spf(:,1),'q',time_spf(:,2));
time_spf = datenum(time_spf,'yyyyqq');
timestr_spf = datestr(time_spf,'yyyy-qq');

[time, ind] = unique(time_spf);
timestr = datestr(time,'yyyy-qq');
T = numel(time);

ids10 = num(start_i:end,3);
spf10 = num(start_i:end,5);


%% Read in individual SR forecasts from SPF, PCE
% (annual rates, %, quarterly frequency)
% Q/Q Rate of Change in the Quarterly- Average Headline CPI Level (annualized percentage points)
xlsx_file = '/Users/lauragati/Dropbox/BC_Research/next/data/raw/SPF/individual_PCE.xlsx';
[num,txt,raw] = xlsread(xlsx_file);
ids1 = num(start_i:end,3);
spf1 = num(start_i:end,10); % CPI1 = yesterday's, CPI2 nowcast, CPI3 one-quarter ahead etc... CPI6 = 4-quarter ahead

%% Read in individual SR forecasts from SPF, PCE core
% (annual rates, %, quarterly frequency)
% Q/Q Rate of Change in the Quarterly- Average Headline CPI Level (annualized percentage points)
xlsx_file = '/Users/lauragati/Dropbox/BC_Research/next/data/raw/SPF/individual_COREPCE.xlsx';
[num,txt,raw] = xlsread(xlsx_file);
ids1 = num(start_i:end,3);
spf1_core = num(start_i:end,10); % CPI1 = yesterday's, CPI2 nowcast, CPI3 one-quarter ahead etc... CPI6 = 4-quarter ahead

%% Need to fish out time series for each ID

ids10_unique = unique(ids10);
ids1_unique = unique(ids1);

% Is it the same guys? 
find(ids10_unique - ids1_unique) % yes
find(ids10 - ids1) % yes totally

% Ok cool then can simplify
ids = ids10;
ids_uni = ids10_unique;
n_fcsters = numel(ids10_unique);

LRE = nan(T,n_fcsters);
SRE = nan(T,n_fcsters);
SRE_core = nan(T,n_fcsters);

n_occur = zeros(n_fcsters,1);
tic
for f95 = 1:n_fcsters
    id = ids_uni(f95);
    idx = find(ids==id);
    n_occur(f95) = numel(idx); 
    
    % now find at what times this forecaster was in the sample
    time_f = time_spf(idx);
    c = ismember(time, time_f);
    time_idx = find(c);
    
    SRE(time_idx,f95) = spf1(idx);
    SRE_core(time_idx,f95) = spf1_core(idx);
    LRE(time_idx,f95) = spf10(idx);

end
toc

% Cut off last periods (b/c infl NaN in the last and we lose a period b/c fe)
infl_pce = infl_pce(1:end-1);
infl_pce_core = infl_pce_core(1:end-1);

SRE = SRE(1:end-2,:);
SRE_core = SRE_core(1:end-2,:);
% LRE = LRE(1:end-2,:);
LRE = LRE(2:end-1,:)-LRE(1:end-2,:);

time = time(1:end-2);
timestr = datestr(time,'yyyy-qq');
T = numel(time);

% Construct SR forecast errors    
fe_pce = infl_pce - SRE;
fe_pce_core = infl_pce_core - SRE_core;

% return

%% Rolling regression with overlapping windows
length_win=20;
n_win = T-length_win;
betahat = zeros(n_win,1);
pvals = zeros(n_win,1);
N  = zeros(n_win,1);
R2 = zeros(n_win,1);
sd = zeros(n_win,1);
ub95 = zeros(n_win,1);
lb95 = zeros(n_win,1);
ub90 = zeros(n_win,1);
lb90 = zeros(n_win,1);

for j=1:T-length_win
    X = vec(fe_pce(j:j+length_win-1, :));
    Y = vec(LRE(j:j+length_win-1, :));
    lm = fitlm(X,Y);
    betahat(j) = lm.Coefficients.Estimate(2);
    pvals(j)   = lm.Coefficients.pValue(2);
    N(j)  = lm.NumObservations;
    R2(j) = lm.Rsquared.Ordinary; % non-adjusted R2
    sd(j) = sqrt(lm.CoefficientCovariance(2,2));
    
%     % upper and lower bounds of 95% CI (assuming beta is normal)
    ub95(j) = betahat(j) +1.96 * sd(j) ;
    lb95(j) = betahat(j) -1.96 * sd(j);
    
    % upper and lower bounds of 90% CI (assuming beta is normal)
    ub90(j) = betahat(j) +1.28 * sd(j) ;
    lb90(j) = betahat(j) -1.28 * sd(j);

end

figspecs = [this_code, '_', nowstr];
[fs, lw] = plot_configs;

time_cropped = linspace(time(1), time(end), n_win);

figure
set(gcf,'color','w'); % sets white background color
set(gcf, 'Position', get(0, 'Screensize')); % sets the figure fullscreen
% plot(time_cropped,ub,'Color', [0.25,0.25,0.25]); hold on
% plot(time_cropped,lb, 'Color', [0.25,0.25,0.25])
h = plot(time_cropped,betahat, 'k','linewidth', lw);
hold on
% x=0:0.01:2*pi;                  %#initialize x array
% y1=sin(x);                      %#create first curve
% y2=sin(x)+.5;                   %#create second curve
X=[time_cropped,fliplr(time_cropped)];                %#create continuous x value array for plotting
Y=[ub95',fliplr(lb95')];              %#create y values for out and then back
f95 = fill(X,Y,[0.25,0.25,0.25],'edgecolor','none');                  %#plot filled area
% Choose a number between 0 (invisible) and 1 (opaque) for facealpha.  
set(f95,'facealpha',.25)

X=[time_cropped,fliplr(time_cropped)];                %#create continuous x value array for plotting
Y=[ub90',fliplr(lb90')];              %#create y values for out and then back
f90 = fill(X,Y,[0.25,0.25,0.25],'edgecolor','none');                  %#plot filled area
% Choose a number between 0 (invisible) and 1 (opaque) for facealpha.  
set(f90,'facealpha',.45)

plot(time_cropped,0*ones(1,length(betahat)), 'k--', 'linewidth',lw)
ax = gca; % current axes
ax.FontSize = fs;
set(gca,'TickLabelInterpreter', 'latex');
grid on
grid minor
legend([h,f95, f90], '$\hat{\beta}_1^w$', '95\% confidence interval','90\% confidence interval','location', 'southoutside', 'interpreter', 'latex', 'NumColumns',3)
legend('boxoff')
datetick('x','yyyy', 'keeplimits')
% The next three lines force the figure to start where the data starts
xaxislimits= get(gca,'XLim');
xaxislimits(1) = time_cropped(1);
xaxislimits(2) = time_cropped(end);
set(gca, 'XLim', xaxislimits);

figname = ['rolling_overlapping_pce_', figspecs];
if print_figs ==1
    disp(figname)
    cd(figpath)
    export_fig(figname)
    cd(current_dir)
    close
end

%% Same thing just regressing on PCE inflation

infl_rep = repmat(infl_pce,1,n_fcsters);
for j=1:T-length_win
%     X = vec(fe(j:j+length_win-1, :));
    X = vec(infl_rep(j:j+length_win-1, :));
    Y = vec(LRE(j:j+length_win-1, :));
    lm = fitlm(X,Y);
    betahat(j) = lm.Coefficients.Estimate(2);
    pvals(j)   = lm.Coefficients.pValue(2);
    N(j)  = lm.NumObservations;
    R2(j) = lm.Rsquared.Ordinary; % non-adjusted R2
    sd(j) = sqrt(lm.CoefficientCovariance(2,2));
    
%     % upper and lower bounds of 95% CI (assuming beta is normal)
    ub95(j) = betahat(j) +1.96 * sd(j) ;
    lb95(j) = betahat(j) -1.96 * sd(j);
    
    % upper and lower bounds of 90% CI (assuming beta is normal)
    ub90(j) = betahat(j) +1.28 * sd(j) ;
    lb90(j) = betahat(j) -1.28 * sd(j);

end

figspecs = [this_code, '_', nowstr];
[fs, lw] = plot_configs;

time_cropped = linspace(time(1), time(end), n_win);

figure
set(gcf,'color','w'); % sets white background color
set(gcf, 'Position', get(0, 'Screensize')); % sets the figure fullscreen
% plot(time_cropped,ub,'Color', [0.25,0.25,0.25]); hold on
% plot(time_cropped,lb, 'Color', [0.25,0.25,0.25])
h = plot(time_cropped,betahat, 'k','linewidth', lw);
hold on
% x=0:0.01:2*pi;                  %#initialize x array
% y1=sin(x);                      %#create first curve
% y2=sin(x)+.5;                   %#create second curve
X=[time_cropped,fliplr(time_cropped)];                %#create continuous x value array for plotting
Y=[ub95',fliplr(lb95')];              %#create y values for out and then back
f95 = fill(X,Y,[0.25,0.25,0.25],'edgecolor','none');                  %#plot filled area
% Choose a number between 0 (invisible) and 1 (opaque) for facealpha.  
set(f95,'facealpha',.25)

X=[time_cropped,fliplr(time_cropped)];                %#create continuous x value array for plotting
Y=[ub90',fliplr(lb90')];              %#create y values for out and then back
f90 = fill(X,Y,[0.25,0.25,0.25],'edgecolor','none');                  %#plot filled area
% Choose a number between 0 (invisible) and 1 (opaque) for facealpha.  
set(f90,'facealpha',.45)

plot(time_cropped,0*ones(1,length(betahat)), 'k--', 'linewidth',lw)
ax = gca; % current axes
ax.FontSize = fs;
set(gca,'TickLabelInterpreter', 'latex');
grid on
grid minor
legend([h,f95, f90], '$\hat{\beta}_1^w$', '95\% confidence interval','90\% confidence interval','location', 'southoutside', 'interpreter', 'latex', 'NumColumns',3)
legend('boxoff')
datetick('x','yyyy', 'keeplimits')
% The next three lines force the figure to start where the data starts
xaxislimits= get(gca,'XLim');
xaxislimits(1) = time_cropped(1);
xaxislimits(2) = time_cropped(end);
set(gca, 'XLim', xaxislimits);

figname = ['rolling_overlapping_pi_pce_', figspecs];
if print_figs ==1
    disp(figname)
    cd(figpath)
    export_fig(figname)
    cd(current_dir)
    close
end

%% PCE CORE

%% Rolling regression with overlapping windows, PCE CORE
length_win=20;
n_win = T-length_win;
betahat = zeros(n_win,1);
pvals = zeros(n_win,1);
N  = zeros(n_win,1);
R2 = zeros(n_win,1);
sd = zeros(n_win,1);
ub95 = zeros(n_win,1);
lb95 = zeros(n_win,1);
ub90 = zeros(n_win,1);
lb90 = zeros(n_win,1);

for j=1:T-length_win
    X = vec(fe_pce_core(j:j+length_win-1, :));
    Y = vec(LRE(j:j+length_win-1, :));
    lm = fitlm(X,Y);
    betahat(j) = lm.Coefficients.Estimate(2);
    pvals(j)   = lm.Coefficients.pValue(2);
    N(j)  = lm.NumObservations;
    R2(j) = lm.Rsquared.Ordinary; % non-adjusted R2
    sd(j) = sqrt(lm.CoefficientCovariance(2,2));
    
%     % upper and lower bounds of 95% CI (assuming beta is normal)
    ub95(j) = betahat(j) +1.96 * sd(j) ;
    lb95(j) = betahat(j) -1.96 * sd(j);
    
    % upper and lower bounds of 90% CI (assuming beta is normal)
    ub90(j) = betahat(j) +1.28 * sd(j) ;
    lb90(j) = betahat(j) -1.28 * sd(j);

end

figspecs = [this_code, '_', nowstr];
[fs, lw] = plot_configs;

time_cropped = linspace(time(1), time(end), n_win);

figure
set(gcf,'color','w'); % sets white background color
set(gcf, 'Position', get(0, 'Screensize')); % sets the figure fullscreen
% plot(time_cropped,ub,'Color', [0.25,0.25,0.25]); hold on
% plot(time_cropped,lb, 'Color', [0.25,0.25,0.25])
h = plot(time_cropped,betahat, 'k','linewidth', lw);
hold on
% x=0:0.01:2*pi;                  %#initialize x array
% y1=sin(x);                      %#create first curve
% y2=sin(x)+.5;                   %#create second curve
X=[time_cropped,fliplr(time_cropped)];                %#create continuous x value array for plotting
Y=[ub95',fliplr(lb95')];              %#create y values for out and then back
f95 = fill(X,Y,[0.25,0.25,0.25],'edgecolor','none');                  %#plot filled area
% Choose a number between 0 (invisible) and 1 (opaque) for facealpha.  
set(f95,'facealpha',.25)

X=[time_cropped,fliplr(time_cropped)];                %#create continuous x value array for plotting
Y=[ub90',fliplr(lb90')];              %#create y values for out and then back
f90 = fill(X,Y,[0.25,0.25,0.25],'edgecolor','none');                  %#plot filled area
% Choose a number between 0 (invisible) and 1 (opaque) for facealpha.  
set(f90,'facealpha',.45)

plot(time_cropped,0*ones(1,length(betahat)), 'k--', 'linewidth',lw)
ax = gca; % current axes
ax.FontSize = fs;
set(gca,'TickLabelInterpreter', 'latex');
grid on
grid minor
legend([h,f95, f90], '$\hat{\beta}_1^w$', '95\% confidence interval','90\% confidence interval','location', 'southoutside', 'interpreter', 'latex', 'NumColumns',3)
legend('boxoff')
datetick('x','yyyy', 'keeplimits')
% The next three lines force the figure to start where the data starts
xaxislimits= get(gca,'XLim');
xaxislimits(1) = time_cropped(1);
xaxislimits(2) = time_cropped(end);
set(gca, 'XLim', xaxislimits);

figname = ['rolling_overlapping_pce_core_', figspecs];
if print_figs ==1
    disp(figname)
    cd(figpath)
    export_fig(figname)
    cd(current_dir)
    close
end

%% Same thing just regressing on PCE CORE inflation

infl_rep = repmat(infl_pce_core,1,n_fcsters);
for j=1:T-length_win
%     X = vec(fe(j:j+length_win-1, :));
    X = vec(infl_rep(j:j+length_win-1, :));
    Y = vec(LRE(j:j+length_win-1, :));
    lm = fitlm(X,Y);
    betahat(j) = lm.Coefficients.Estimate(2);
    pvals(j)   = lm.Coefficients.pValue(2);
    N(j)  = lm.NumObservations;
    R2(j) = lm.Rsquared.Ordinary; % non-adjusted R2
    sd(j) = sqrt(lm.CoefficientCovariance(2,2));
    
%     % upper and lower bounds of 95% CI (assuming beta is normal)
    ub95(j) = betahat(j) +1.96 * sd(j) ;
    lb95(j) = betahat(j) -1.96 * sd(j);
    
    % upper and lower bounds of 90% CI (assuming beta is normal)
    ub90(j) = betahat(j) +1.28 * sd(j) ;
    lb90(j) = betahat(j) -1.28 * sd(j);

end

figspecs = [this_code, '_', nowstr];
[fs, lw] = plot_configs;

time_cropped = linspace(time(1), time(end), n_win);

figure
set(gcf,'color','w'); % sets white background color
set(gcf, 'Position', get(0, 'Screensize')); % sets the figure fullscreen
% plot(time_cropped,ub,'Color', [0.25,0.25,0.25]); hold on
% plot(time_cropped,lb, 'Color', [0.25,0.25,0.25])
h = plot(time_cropped,betahat, 'k','linewidth', lw);
hold on
% x=0:0.01:2*pi;                  %#initialize x array
% y1=sin(x);                      %#create first curve
% y2=sin(x)+.5;                   %#create second curve
X=[time_cropped,fliplr(time_cropped)];                %#create continuous x value array for plotting
Y=[ub95',fliplr(lb95')];              %#create y values for out and then back
f95 = fill(X,Y,[0.25,0.25,0.25],'edgecolor','none');                  %#plot filled area
% Choose a number between 0 (invisible) and 1 (opaque) for facealpha.  
set(f95,'facealpha',.25)

X=[time_cropped,fliplr(time_cropped)];                %#create continuous x value array for plotting
Y=[ub90',fliplr(lb90')];              %#create y values for out and then back
f90 = fill(X,Y,[0.25,0.25,0.25],'edgecolor','none');                  %#plot filled area
% Choose a number between 0 (invisible) and 1 (opaque) for facealpha.  
set(f90,'facealpha',.45)

plot(time_cropped,0*ones(1,length(betahat)), 'k--', 'linewidth',lw)
ax = gca; % current axes
ax.FontSize = fs;
set(gca,'TickLabelInterpreter', 'latex');
grid on
grid minor
legend([h,f95, f90], '$\hat{\beta}_1^w$', '95\% confidence interval','90\% confidence interval','location', 'southoutside', 'interpreter', 'latex', 'NumColumns',3)
legend('boxoff')
datetick('x','yyyy', 'keeplimits')
% The next three lines force the figure to start where the data starts
xaxislimits= get(gca,'XLim');
xaxislimits(1) = time_cropped(1);
xaxislimits(2) = time_cropped(end);
set(gca, 'XLim', xaxislimits);

figname = ['rolling_overlapping_pi_pce_core', figspecs];
if print_figs ==1
    disp(figname)
    cd(figpath)
    export_fig(figname)
    cd(current_dir)
    close
end


    