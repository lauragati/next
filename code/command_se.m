function se = command_se
% command_se.m
% Compute standard errors of bootstrapped alpha_hat
% 6 June 2021

clearvars
close all
clc

datestr(now)

%% Load estimated alpha_hat, should be 5x100 (nfe x nboot)
% filename = alpha_hat_boot_date.mat;
% load([filename, '.mat'])

%% Compute standard errors

% Pretend like you had the estimates:
rng(0)
alpha_hat = rand(5,100);
[~,NB] = size(alpha_hat);

% Compute std.dev with Matlab's inbuilt function
stddev_matlab = std(alpha_hat,0,2);

% And your own:
stddev_me = sqrt(1/(NB-1) .* sum((mean(alpha_hat,2) - alpha_hat).^2,2)  );

disp('If these aren''t zeros, there''s a problem:')
stddev_matlab - stddev_me % ok cool they're equal.


% Then the std error is the stddev/sqrt(n)
se = stddev_me/sqrt(NB);

% Now you just need to add these to the plot as alpha_hat +- se