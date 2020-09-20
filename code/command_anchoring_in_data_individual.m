% command_anchoring_in_data_individual.m
% Look at evidence of unanchoring in data from individual forecasters (following Jenny Tang)
% 20 Sept 2020

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

%% Read in CPI
% observation_start = '1991-10-01'; % LR-E is the binding constraint at the beginning.
observation_start = '1991-01-01'; % need to take a little earlier to construct inflation and forecast errors
observation_end   = '2020-07-01'; 
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

%% Read in individual LR forecasts from SPF
% (annual rates, %, quarterly frequency)
% Ten-Year Annual-Average Headline CPI Inflation (Horizon: Year of Survey and Following Nine Years)
xlsx_file = '/Users/lauragati/Dropbox/BC_Research/next/data/raw/SPF/individual_CPI10.xlsx';
[num,txt,raw] = xlsread(xlsx_file);
start_i = 3814-1;
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


%% Read in individual SR forecasts from SPF
% (annual rates, %, quarterly frequency)
% Q/Q Rate of Change in the Quarterly- Average Headline CPI Level (annualized percentage points)
xlsx_file = '/Users/lauragati/Dropbox/BC_Research/next/data/raw/SPF/individual_CPI.xlsx';
[num,txt,raw] = xlsread(xlsx_file);
ids1 = num(start_i:end,3);
spf1 = num(start_i:end,10); % CPI1 = yesterday's, CPI2 nowcast, CPI3 one-quarter ahead etc... CPI6 = 4-quarter ahead

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
n_occur = zeros(n_fcsters,1);
tic
for f = 1:n_fcsters
    id = ids_uni(f);
    idx = find(ids==id);
    n_occur(f) = numel(idx); 
    
    % now find at what times this forecaster was in the sample
    time_f = time_spf(idx);
    c = ismember(time, time_f);
    time_idx = find(c);
    
    SRE(time_idx,f) = spf1(idx);
    LRE(time_idx,f) = spf10(idx);

end
toc

% Cut off last periods (b/c infl NaN in the last and we lose a period b/c fe)
infl = infl(1:end-1);
SRE = SRE(1:end-2,:);
% LRE = LRE(1:end-2,:);
LRE = LRE(2:end-1,:)-LRE(1:end-2,:);

time = time(1:end-2);
timestr = datestr(time,'yyyy-qq');
T = numel(time);

% Construct SR forecast errors    
fe = infl - SRE;

%% One big regression

fitlm(vec(fe),vec(LRE))


%% % Let's implement a rolling regression
n_windows = 6; % 3,6,19, 38
length_window = T/n_windows;
idx_window = 1:T/n_windows:T;
windows = datestr(time(idx_window),'yyyy-qq');
% add the last
idx_window = [idx_window, T];

bet_win = zeros(n_windows,1);
betahat = zeros(n_windows,1);
pvals = zeros(n_windows,1);
N  = zeros(n_windows,1);
R2 = zeros(n_windows,1);

for i=1:n_windows
    X = vec(fe(idx_window(i)+1:idx_window(i+1), :));
    Y = vec(LRE(idx_window(i)+1:idx_window(i+1),:));
%     X = abs(vec(fe(idx_window(i)+1:idx_window(i+1), :)));
%     Y = abs(vec(LRE(idx_window(i)+1:idx_window(i+1),:)));
    bet_win(i) = (X'*X) \ (X'*Y); % these become nans cause matlab gets screwed up with nan-multiplication
    lm = fitlm(X,Y);
    betahat(i) = lm.Coefficients.Estimate(2);
    pvals(i)   = lm.Coefficients.pValue(2);
    N(i)  = lm.NumObservations;
    R2(i) = lm.Rsquared.Ordinary; % non-adjusted R2
end

windows(1:n_windows, :)
[betahat, pvals, R2]
