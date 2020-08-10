% commmand_verify_VAR_OLS.m
% check that the VAR estimation codes are correct
% 5 August 2020

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

%% Compute weighting matrix and initialize alpha

filename = 'acf_sim_univariate_data_06_Jul_2020'; % simulated data, nfe=5, fe=(-2,2), alph_true = (0.05; 0.025; 0; 0.025; 0.05); see Notes 6 July 2020

load([filename, '.mat'])
nobs = acf_outputs{3};
p = acf_outputs{4};
filt_data = acf_outputs{6};

% First make sure that having a nvar-dimensional Y is the same as
% estimating equation by equation, and it is
[~,B,res,sigma] = sr_var(filt_data', p);
[B_RFee,res_RFee,sigma_RFee] = rf_var_eq_by_eq(filt_data', p);
B - B_RFee % same
res - res_RFee % same
sigma - sigma_RFee % same

% Now make sure that the regressor is computed correctly
tic
[B_RF,res_RF,sigma_RF] = reduform_var(filt_data', p);
toc
tic
[B_newreg,res_newreg,sigma_newreg] = rf_var(filt_data', p); % mine is slower unfortunately but at least it avoids eval
toc
find(B_RF - B_newreg) % same
find(res_RF - res_newreg) % same
find(sigma_RF - sigma_newreg) % same

% But the point is that the old codes, sr_var and reduform_var are now
% verified to be working correctly. 

