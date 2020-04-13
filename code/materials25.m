% materials25
% Prepare macro lunch of 15 April 2020
% --> get some kind of implementation of Taylor rule to work
% 1. Return to fsolve with the cleaned up learning-simulation code
% that doesn't work

clearvars
close all
clc

% Add all the relevant paths and grab the codename
this_code = mfilename;
[current_dir, basepath, BC_researchpath,toolpath,export_figpath,figpath,tablepath,datapath] = add_paths;
todays_date = strrep(datestr(today), '-','_');

% Variable stuff ---
print_figs        = 0;
stop_before_plots = 0;
skip_old_plots    = 0;
output_table = print_figs;

skip = 1;

[fs, lw, grey, silver, maroon, light_coral, light_salmon,dark_green, green, light_green, light_sky_blue, ...
    teal, purple, saddle_brown, light_brown, fs_pres,lw_pres] = plot_configs;

%% do an inflation expectation plot for prezi
if skip==0
% UMich Inflation Expectations
series_id = 'MICH';
observation_start = '1978-01-01';
observation_end   = '1990-01-01';
[output] = getFredData(series_id, observation_start, observation_end);
umich = output.Data(:,2);
time_umich = output.Data(:,1);


figure
set(gcf,'color','w'); % sets white background color
set(gcf, 'Position', get(0, 'Screensize')); % sets the figure fullscreen
plot(time_umich,umich, 'k', 'linewidth',lw_pres)
ax = gca; % current axes
ax.FontSize = fs_pres;
grid on
grid minor
datetick('x','yyyy')
if print_figs ==1
    figname = ['macrolunch_umich']
    cd(figpath)
    export_fig(figname)
    cd(current_dir)
    close
end

end

%% 

T = 100
N = 10
ndrop =0; ne=3;
[param, set, param_names, param_values_str, param_titles] = parameters_next;

sig_r = param.sig_r;
sig_i = param.sig_i;
sig_u = param.sig_u;

% Params for the general learning code
constant_only = 1; % learning constant only (vector - _smooth won't work in this case)
constant_only_pi_only = 11; % learning constant only, inflation only (scalar)
slope_and_constant = 2;

dgain = 1;  % 1 = decreasing gain, 2 = endogenous gain, 3 = constant gain
again_critCEMP  = 21;
again_critCUSUM = 22;
again_critsmooth = 23;
cgain = 3;

% Model selection
%%%%%%%%%%%%%%%%%%%
PLM = constant_only;
gain = again_critsmooth;
%%%%%%%%%%%%%%%%%%%

% RE model
[fyn, fxn, fypn, fxpn] = model_NK(param);
[gx,hx]=gx_hx_alt(fyn,fxn,fypn,fxpn);
[ny, nx] = size(gx);
[Aa, Ab, As] = matrices_A_13_true_baseline(param, hx);
SIG = eye(nx).*[sig_r, sig_i, sig_u]';
eta = SIG; %just so you know

rng(0)
eN = randn(ne,T+ndrop,N);
e = squeeze(eN(:,:,1));
eN(2,:,:) = zeros(T+ndrop,N);

% an initial simulation
[PLM_name, gain_name, gain_title] = give_names(PLM, gain);
[xsim, ysim, k, phi_seq, FA, FB, diff] = sim_learnLH_clean_g(param,gx,hx,eta, Aa, Ab, As,PLM, T,ndrop,e);

% create_plot_observables(ysim)
% create_plot_observables(1./k)

% Optimization Parameters
options = optimoptions('fmincon', 'TolFun', 1e-12, 'display', 'iter', 'MaxFunEvals', 10000);

coeffs0 = [param.psi_pi, param.psi_x,param.psi_k,param.psi_pibar, param.psi_xbar];



ub = [5,0,0.9,0.5,0];
% ub = [5,2,0,0,0];% a standard TR, yields coeffs_opt =  (  1.4912    0.0097         0         0         0) for a loss of   1.5981e+03
% ub = [5,0,0,0,0];% a simple TR, yields coeffs_opt =  (1.5287         0         0         0         0) for a loss of   3.8385e+03
ub = [5,0,0.9,0.5,0]; % add k and psibar -> coeffs_opt  (1.5639         0    0.0133    0.0000         0) for a loss of  5.4353e+03
% if you start this same one at lb, you get coeffs_opt  (1        0    0  0.0000         0) for a loss of  138.2685
lb = [1,0,0,0,0];

% let's start at other places
coeffs0 = lb;
loss = objective_ramsey_materials25(coeffs0,param,gx,hx,eta,ndrop,eN)

tic
objh = @(coeffs) objective_ramsey_materials25(coeffs,param,gx,hx,eta,ndrop,eN);
[coeffs_opt, loss_opt] = fmincon(objh, coeffs0, [],[],[],[],lb,ub,[],options)
% fmincon(FUN,X0,A,B,Aeq,Beq,LB,UB,NONLCON,OPTIONS)
toc

param.psi_pi = coeffs_opt(1);
param.psi_x = coeffs_opt(2);
param.psi_k     = coeffs_opt(3);
param.psi_pibar = coeffs_opt(4);
param.psi_xbar  = coeffs_opt(5);

% simulate again
[xsim, ysim, k, phi_seq, FA, FB, diff] = sim_learnLH_clean_g(param,gx,hx,eta, Aa, Ab, As,PLM, T,ndrop,e);
create_plot_observables(ysim)
create_plot_observables(1./k)