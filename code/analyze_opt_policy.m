% analyze_opt_policy.m
% Load value function iteration and some parameterized expectations results
% to analyze the policy function
% 30 May 2020

clearvars
close all
clc

% Add all the relevant paths and grab the codename
this_code = mfilename;
[current_dir, basepath, BC_researchpath,toolpath,export_figpath,figpath,tablepath,datapath, tryouts_path] = add_paths;
todays_date = strrep(datestr(today), '-','_');

% Variable stuff ---
print_figs        = 0;
stop_before_plots = 0;
skip_old_plots    = 0;
output_table = print_figs;

skip = 1;

%% Comparative statics of policy wrt k1 and pibar
[param, set, param_names, param_values_str, param_titles] = parameters_next;
sig_r = param.sig_r;
sig_u = param.sig_u;

% Grids
ns = 2;
sgrid = linspace(-sig_r,sig_r,ns);

% value_output_name = 'value_outputs_approx30_Jun_2020_10_19_45'; % univariate approximated LOM gain
% value_output_name = 'value_outputs_approx17_Jul_2020_16_37_34'; % univariate approximated LOM gain % materials 37
% value_output_name = 'value_outputs_approx27_Aug_2020_14_18_38'; % Calibration C of Materials 43;  pgrid = linspace(-0.1,0.1,np);
% value_output_name = 'value_outputs_approx27_Aug_2020_14_28_32'; % Calibration C of Materials 43;  pgrid = linspace(-0.2,0.2,np); % nicer int-rate magnitudes
% value_output_name = 'value_outputs_approx12_Sep_2020_09_29_49'; % Sept 15 estimation of Materials 44; pgrid = linspace(-1,1,np); % honestly I think they're exactly the same
value_output_name = 'value_outputs_approx17_Sep_2020_14_01_16'; % Sept 21 draft Materials 44;  pgrid = linspace(-0.2,0.2,np);

load([value_output_name, '.mat'])

pp     = value_sols{1};
v      = value_sols{2};
it     = value_sols{3};
pibp   = value_sols{4};
k1p    = value_sols{5};
pgrid  = value_sols{6};
% k1grid = value_sols{7};

% take the history of states from parametric expectations
% pea_output_name = 'pea_outputs_approx30_Jun_2020_09_39_12'; % univariate approximated LOM gain
% pea_output_name = 'pea_outputs_approx17_Jul_2020_11_40_09'; % univariate approximated LOM gain, improved estimation, rng(0) % materials 37
% pea_output_name = 'pea_outputs_approx23_Aug_2020_14_38_14'; % 23 August 2020 calibration
% pea_output_name = 'pea_outputs_approx27_Aug_2020_14_45_53';  % Calibration C of Materials 43; rng(2) default
% pea_output_name = 'pea_outputs_approx27_Aug_2020_14_56_40';  % Calibration C of Materials 43; rng(3)
% pea_output_name = 'pea_outputs_approx27_Aug_2020_15_00_03';  % Calibration C of Materials 43; rng(4)
% pea_output_name = 'pea_outputs_approx12_Sep_2020_09_15_24';  % Sept 15 estimation of Materials 44; rng(2) default
pea_output_name = 'pea_outputs_approx17_Sep_2020_13_47_33'; % Sept 21 draft Materials 44; rng(2)






% load('inputs.mat')
load([pea_output_name, '.mat'])
e = output{1};
ysim7 = output{2};
k7    = output{3};
phi7  = output{4};
seq_opt = output{5};
i_pe = seq_opt(3,:);

% compile state vector
T=length(e)-2; % drop the first and last obs
k1sim = 1./k7(2:end-1);
pibsim = squeeze(phi7(1,1,2:end-1))';

% cheeck max, mean and min in simulation
% max(k1sim) % 0.2699
% mean(k1sim) % 0.0092
% min(k1sim) % 7.7349e-04
% max(pibsim) % 0.0156, 0.2136
% mean(pibsim) % -0.0012, 0.1945
% min(pibsim) % -0.0438, -0.0931

% get 10th and 90th percentile of k1 (measures of anchored and unanchored)
k1sim_sort = sort(k1sim);
if mod(T,10)~=0
    k1_10 = k1sim_sort(floor(T/10)+1);
    k1_90 = k1sim_sort(ceil(9*T/10)+1);
elseif mod(T,10)==0
    k1_10 = (k1sim_sort(floor(T/10)) + k1sim_sort(floor(T/10)+1)) /2;
    k1_90 = (k1sim_sort(ceil(9*T/10)) + k1sim_sort(ceil(9*T/10)+1)) /2;
end



ssim = e(:,2:end-1);


n = 10;
pibvals = linspace(-0.1, 0.1, n); % most values in simulations are between these
pmean = zeros(1,n);
smean = zeros(1,n);

% approximate policy function
ppi = csapi({pgrid,sgrid,sgrid,sgrid,sgrid},it);

%hold every state at their means, just move pibar
i_p = fnval(ppi,[pibvals;smean;smean;smean;smean]);

% % Annualization of inflation and interest rate - I won't do it b/c it
% would need to scale up pibvals too
% i_p =  ((i_p/100+1).^4 -1)*100;
% pibvals =  ((pibvals/100+1).^4 -1)*100;

xlplus = [0.1, 0.3];
ylplus = [0.02, 6.1]; 
create_pretty_plot_x(pibvals,i_p,'$\bar{\pi}$','$i(\bar{\pi}, \cdot)$',xlplus, ylplus,[this_code, '_ip', todays_date],print_figs)

return

%% Regress i on pi and x - this is just to see that under optimal policy, the CB responds much more to inflation and x than under TR
X = [ones(size(i_pe)); ysim7(1,2:end-1); ysim7(2,2:end-1)]';
Y = i_pe';
bet_ols = (X'*X)^(-1)*(X'*Y);
yhat = X*bet_ols;
ybar = mean(Y);
% sum of squared residuals:
SSR = sum((Y-yhat).^2);
% total sum of squares:
SST = sum((Y-ybar).^2);
R2 = 1 - SSR/SST;

mdl = fitlm(X(:,2:end),Y)


%% Regress pibar on k1, pibar{t-1}, fe and shocks - to see how much movement in LR exp is associated with what gains
pibt   = pibsim(2:end);
pibt_1 = pibsim(1:end-1);
k1t    = k1sim(2:end);
fe     = ysim7(1,3:end-1) - squeeze(phi7(1,1,3:end-1))';
shocks = e([1;3],3:end-1);

X = [ones(size(k1t)); k1t; pibt_1; fe]';
Y = pibt';
bet_ols = (X'*X)^(-1)*(X'*Y);
yhat = X*bet_ols;
ybar = mean(Y);
% sum of squared residuals:
SSR = sum((Y-yhat).^2);
% total sum of squares:
SST = sum((Y-ybar).^2);
R2 = 1 - SSR/SST;

mdl = fitlm(X(:,2:end),Y)

% Calculation 30 July 2020
% So if pibt_1 = 0, and forecast error increases by 1,
% we know from the estimated figure that the gain increases to 0.02, so the
% change in pibar wrt fe is:
betPibFe = 0.02*bet_ols(2)+bet_ols(4);


