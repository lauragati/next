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
XX = X(:,2:end); % without constant

mdl = fitlm(XX,Y)

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

B1 = ridge(Y,XX,lam); % Matlab uses standardized X and demeaned Y for any choice of scaled, 0 or 1.
B0 = ridge(Y,XX,lam,0); % When scaled is 1 (default), ridge does not restore the coefficients to the original data scale. When scaled is 0, ridge restores the coefficients to the scale of the original data. 
bet_ridge = (X'*X+ lam*eye(p))^(-1)*(X'*Y);

results_mat1 = {'OLS', 'Ridge', 'Ridge Matlab, scale=0', 'Ridge Matlab, scale=1'} 
[bet_ols, bet_ridge, B0, [0;B1]]
disp('Ridge Matlab, scale=0 seems the way to go')
disp('I think that ''Ridge Matlab, scale=1'' implies a zero intercept by assumption b/c X has mean 0')

%% Ridge when setting lambda>0 

lam = 0.1;

B1 = ridge(Y,XX,lam);
B0 = ridge(Y,XX,lam,0);
bet_ridge = (X'*X+ lam*eye(p))^(-1)*(X'*Y);

results_mat2 = {'OLS', 'Ridge', 'Ridge Matlab, scale=0', 'Ridge Matlab, scale=1'} 
[bet_ols, bet_ridge, B0, [0;B1]]
disp('I should be getting the same thing as ''Ridge Matlab, scale=1'' ')

%% Try ridge again, this time standardizing (demean Y, standardize X)

% Yhat is the centered (demeaned) Y
Yhat = (Y-mean(Y));

% denote by hat the standardized regressor
Xhat = (X-mean(X))./std(X);
% Take out intercept, it's 0/0 anyway
Xhat = Xhat(:,2:end);
p = size(Xhat,2);

B1 = ridge(Y,XX,lam);
B0 = ridge(Y,XX,lam,0);
% Note: Matlab doesn't use the standardized Y
bet_ridge = (Xhat'*Xhat+ lam*eye(p))^(-1)*(Xhat'*Y);

results_mat3 = {'OLS', 'Ridge', 'Ridge Matlab, scale=0', 'Ridge Matlab, scale=1'} 
[bet_ols, [0;bet_ridge], B0, [0;B1]]

% scale B1 back as Matlab does
B0 - [mean(Y)-mean(XX)*(B1./std(XX,0,1)'); B1./std(XX,0,1)']


%% Investigate how to scale back depending on whether you demean Y or not

B1_w = (Xhat'*Xhat+ lam*eye(p))^(-1)*(Xhat'*Yhat);
B1_wo = (Xhat'*Xhat+ lam*eye(p))^(-1)*(Xhat'*Y);
% you get the same estimate -> See notes 18 August 2020

B_rescaled_w = [mean(Y)-mean(XX)*(B1_w./std(XX,0,1)'); B1_w./std(XX,0,1)'];

% I would have thought that if you don't use demeaned Y, then when scaling back, you don't put mean(Y) back in the intercept:
B_rescaled_wo = [-mean(XX)*(B1_wo./std(XX,0,1)'); B1_wo./std(XX,0,1)'];
% But it turns out that's wrong. B/c B1_wo = B1_w, you don't get the same.

% Conclusion: add mean(Y) back when rescaling in any case. It doesn't
% matter if you demean Y or not, so easiest may be not to demean.

%% Multidimensional Y - ridge for VAR

nlags = 2;
dataset = [gini, gdp, u];

T = size(dataset,1); % time periods
nvar = size(dataset,2);

% Create a dataset with lags of all variables
% time, variable, lag
lagged_data = zeros(T-nlags,nvar,nlags+1);
for i=1:nlags+1
    p=i-1;
    lag_p = dataset(nlags+1-p:end-p,:);
    lagged_data(:,:,i) = lag_p;
end

% Regressand
Y = squeeze(lagged_data(:,:,1));
% Regressor
X = reshape(lagged_data(:,:,2:end),T-nlags,nvar*nlags);
% Don't add constant for ridge

% standardize X
Xhat = (X-mean(X))./std(X);
p = size(Xhat,2);
bet_ridge = (Xhat'*Xhat+ lam*eye(p))^(-1)*(Xhat'*Y);
% scale back 
bet_ridge_rescaled = [mean(Y)-mean(X)*(bet_ridge./std(X,0,1)'); bet_ridge./std(X,0,1)'];

% Now let's see if Matlab can do it
% B0 = ridge(Y,X,lam,0); % "Y must be a column vector" error
B0 = zeros(size(bet_ridge_rescaled));
for i=1:3
    B0(:,i) = ridge(Y(:,i),X,lam,0);
end

% the grand finale: you get the same to e-09
bet_ridge_rescaled - B0
    