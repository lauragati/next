% compare_value_pea_results_approx.m
% Load value function iteration and parameterized expectations results for the approximated version of the LOM gain and plot them
% to compare
% 13 June 2020

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
stop_before_plots = 0;
skip_old_plots    = 0;
output_table = print_figs;

skip = 1;

%%
[param, set, param_names, param_values_str, param_titles] = parameters_next;
sig_r = param.sig_r;
sig_u = param.sig_u;

% RE model
[fyn, fxn, fypn, fxpn] = model_NK(param);
[gx,hx]=gx_hx_alt(fyn,fxn,fypn,fxpn);

% Grids
% nk = 4;
% gbar = param.gbar;
% k1grid = linspace(0,gbar,nk);
% np = 4;
% pgrid = linspace(-0.2,0.2,np);
% pgrid = linspace(-10,10,np); % this fucker causes the difference!
% pgrid = linspace(-4,4,np);


ns = 2;
sgrid = linspace(-sig_r,sig_r,ns);
p = 0.5;
PI = [p*p, p*(1-p); (1-p)*p, (1-p)*(1-p)];

% value_output_name = 'value_outputs_approx30_Jun_2020_10_19_45'; % univariate approximated LOM gain
% value_output_name = 'value_outputs_approx17_Jul_2020_14_27_00'; % univariate approximated LOM gain, improved estimation, pbar in (-1,1), np=4
% value_output_name = 'value_outputs_approx17_Jul_2020_16_13_24'; % univariate approximated LOM gain, improved estimation, pbar in (-1,1), np=6
% value_output_name = 'value_outputs_approx17_Jul_2020_16_37_34'; % univariate approximated LOM gain, improved estimation, pbar in (-1,1), np=8
% value_output_name = 'value_outputs_approx23_Aug_2020_14_50_26'; % 23 August 2020 calibration 
% value_output_name = 'value_outputs_approx27_Aug_2020_14_18_38'; % Calibration C of Materials 43;  pgrid = linspace(-0.1,0.1,np); 
% value_output_name = 'value_outputs_approx27_Aug_2020_14_28_32'; % Calibration C of Materials 43;  pgrid = linspace(-0.2,0.2,np); % nicer int-rate magnitudes
% value_output_name = 'value_outputs_approx12_Sep_2020_09_41_01'; % Sept 15 estimation of Materials 44; pgrid = linspace(-0.1,0.1,np);
% value_output_name = 'value_outputs_approx12_Sep_2020_09_19_47'; % Sept 15 estimation of Materials 44; pgrid = linspace(-0.2,0.2,np);
% value_output_name = 'value_outputs_approx12_Sep_2020_09_29_49'; % Sept 15 estimation of Materials 44; pgrid = linspace(-1,1,np); % honestly I think they're exactly the same
% value_output_name = 'value_outputs_approx17_Sep_2020_13_51_53'; % Sept 21 draft Materials 44;  pgrid = linspace(-0.1,0.1,np);
value_output_name = 'value_outputs_approx17_Sep_2020_14_01_16'; % Sept 21 draft Materials 44;  pgrid = linspace(-0.2,0.2,np);



load([value_output_name, '.mat'])

pp     = value_sols{1};
v      = value_sols{2};
it     = value_sols{3};
pibp   = value_sols{4};
k1p    = value_sols{5};
pgrid  = value_sols{6};
% k1grid = value_sols{7};


%%
% take the history of states from parametric expectations
% pea_output_name = 'pea_outputs_approx30_Jun_2020_09_39_12'; % univariate approximated LOM gain
% pea_output_name = 'pea_outputs_approx17_Jul_2020_11_40_09'; % univariate approximated LOM gain, improved estimation, rng(0)
% pea_output_name = 'pea_outputs_approx17_Jul_2020_15_28_03'; % univariate approximated LOM gain, improved estimation, rng(2)
% pea_output_name = 'pea_outputs_approx23_Aug_2020_14_38_14'; % 23 August 2020 calibration
% pea_output_name = 'pea_outputs_approx27_Aug_2020_14_45_53';  % Calibration C of Materials 43; rng(2) default
% pea_output_name = 'pea_outputs_approx27_Aug_2020_14_56_40';  % Calibration C of Materials 43; rng(3)
% pea_output_name = 'pea_outputs_approx27_Aug_2020_15_00_03';  % Calibration C of Materials 43; rng(4)
% pea_output_name = 'pea_outputs_approx12_Sep_2020_09_15_24';  % Sept 15 estimation of Materials 44; rng(2) default
pea_output_name = 'pea_outputs_approx17_Sep_2020_13_47_33'; % Sept 21 draft Materials 44; rng(2)




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
ssim = e(:,2:end-1);
ssimt_1 = e(:,1:end-2);
X = [pibsim; ssim(1,:); ssim(3,:); ssimt_1(1,:); ssimt_1(3,:)];

% approximate policy function
% ppi = csapi({k1grid,pgrid,sgrid,sgrid,sgrid,sgrid},it);
ppi = csapi({pgrid,sgrid,sgrid,sgrid,sgrid},it);
i_vi = fnval(ppi,X);

policies = [i_pe; i_vi];

% Annualization of inflation and interest rate
ysim7([1,3],:) = ((ysim7([1,3],:)/100+1).^4 -1)*100;
policies =  ((policies/100+1).^4 -1)*100;


% % save some nice plots for draft and prezis
xlab = 'Quarters';
ylab = 'dont display';

xlplus = [45,0.01];
ylplus = [-400,0];
create_pretty_plot_holdon(policies, {'PEA', 'VFI'},xlab,ylab,xlplus, ylplus, [this_code, '_', value_output_name, '_', pea_output_name, '_pretty_', nowstr], print_figs)
create_pretty_subplots(ysim7(:,2:end-1),{'$\pi$', '$x$','$i$'},xlab,ylab,xlplus, ylplus, ['implement_anchTC_obs_approx', nowstr], print_figs)

xlplus = [20,0.01];
ylplus = [-400,0];
create_pretty_plot_holdon(policies(:,1:40), {'PEA', 'VFI'},xlab,ylab,xlplus, ylplus, [this_code, '_', value_output_name, '_', pea_output_name, '_pretty_40q_', nowstr], print_figs)
create_pretty_subplots(ysim7(:,2:40),{'$\pi$', '$x$','$i$'},xlab,ylab,xlplus, ylplus, ['implement_anchTC_obs_approx_40q', nowstr], print_figs)


return
xlplus = [45,0.0025];
ylplus = [-100,0];
create_pretty_plot_x(1:length(k7)-1,1./k7(1:end-1), xlab,ylab,xlplus, ylplus,['implement_anchTC_invgain', nowstr], print_figs)
create_pretty_plot_x(1:T,pibsim,xlab,ylab,xlplus, ylplus,['implement_anchTC_pibar', nowstr], print_figs)


