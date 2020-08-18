% command_verify_ridge.m
% Check on a random (not VAR) estimation that you're doing ridge correctly.
% 17 August 2020


clearvars
close all
clc

% Add all the relevant paths and grab the codename
this_code = mfilename;
[current_dir, basepath, BC_researchpath,toolpath,export_figpath,figpath,tablepath,datapath] = add_paths;
todays_date = strrep(datestr(today), '-','_');
nowstr = strrep(strrep(strrep(datestr(now), '-','_'), ' ', '_'), ':', '_');

% Variable stuff ---
print_figs        = 1;
if contains(current_dir, 'gsfs0') % sirius server
    print_figs=1;
end
stop_before_plots = 0;
skip_old_plots    = 0;
output_table = print_figs;

save_estim_outputs =0;

skip = 1;
investigate_loss=0;
[fs, lw] = plot_configs;
datestr(now)

%% Load data

% Income Gini coefficient for US households
[output8] = getFredData('GINIALLRH'); %, observation_start, observation_end, units, frequency, aggregation_method);
gini = output8.Data(:,2);

range = char(output8.DataRange);
obs_start = range(1:10);
obs_end = range(15:end);

% Real GDP
[output4] = getFredData('GDPC1', obs_start, obs_end, 'lin', 'a', 'avg');
gdp = output4.Data(:,2);


% unemployment
[output6] = getFredData('UNRATE', obs_start, obs_end, 'lin', 'a', 'avg');
u = output6.Data(:,2);

Y = gini;
X = [ones(size(gdp)), gdp, u];

[n,p] = size(X);

%% OLS

mdl = fitlm(X(:,2:end),Y)

bet_ols = (X'*X)^(-1)*(X'*Y);
yhat = X*bet_ols;
ybar = mean(Y);
% sum of squared residuals:
SSR = sum((Y-yhat).^2);
% total sum of squares:
SST = sum((Y-ybar).^2);
R2 = 1 - SSR/SST;

%% Ridge when setting lambda=0 (should give OLS)

lam = 0;

B1 = ridge(Y,X(:,2:end),lam); % Matlab uses standardized X and Y for any choice of scaled, 0 or 1.
B0 = ridge(Y,X(:,2:end),lam,0); % When scaled is 0, ridge restores the coefficients to the scale of the original data. 
bet_ridge = (X'*X+ lam*eye(p))^(-1)*(X'*Y);

results_mat1 = {'OLS', 'Ridge', 'Ridge Matlab, scale=0', 'Ridge Matlab, scale=1'} 
[bet_ols, bet_ridge, B0, [0;B1]]
disp('Ridge Matlab, scale=0 seems the way to go')
disp('I think that ''Ridge Matlab, scale=1'' implies a zero intercept by assumption b/c X has mean 0')

%% Ridge when setting lambda>0 

lam = 0.1;

B1 = ridge(Y,X(:,2:end),lam);
B0 = ridge(Y,X(:,2:end),lam,0);
bet_ridge = (X'*X+ lam*eye(p))^(-1)*(X'*Y);

results_mat2 = {'OLS', 'Ridge', 'Ridge Matlab, scale=0', 'Ridge Matlab, scale=1'} 
[bet_ols, bet_ridge, B0, [0;B1]]
disp('I should be getting the same thing as ''Ridge Matlab, scale=1'' ')

%% Try ridge again, this time standardizing
% denote by hat the standardized variables
sdY = sqrt(var(Y));
Yhat = (Y-mean(Y))/sdY;

sdX = sqrt(var(X)); % checked = std(X)
Xhat = (X-mean(X))./sdX;
% Take out intercept, it's 0/0 anyway
Xhat = Xhat(:,2:end);

B1 = ridge(Y,X(:,2:end),lam);
B0 = ridge(Y,X(:,2:end),lam,0);
% Note: Matlab doesn't use the standardized Y
bet_ridge = (Xhat'*Xhat+ lam*eye(p-1))^(-1)*(Xhat'*Y);

results_mat3 = {'OLS', 'Ridge', 'Ridge Matlab, scale=0', 'Ridge Matlab, scale=1'} 
[bet_ols, [0;bet_ridge], B0, [0;B1]]

