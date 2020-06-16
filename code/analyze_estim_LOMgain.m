% analyze_estim_LOMgain.m
% Load value function iteration and some parameterized expectations results
% to analyze the policy function
% 14 June 2020

clearvars
close all
clc

% Add all the relevant paths and grab the codename
this_code = mfilename;
[current_dir, basepath, BC_researchpath,toolpath,export_figpath,figpath,tablepath,datapath, inputsRyan_path] = add_paths;
todays_date = strrep(datestr(today), '-','_');

% Variable stuff ---
print_figs        = 0;
stop_before_plots = 0;
skip_old_plots    = 0;
output_table = print_figs;

skip = 1;
[fs, lw] = plot_configs;


% Load estimation outputs
% filename = 'estim_LOMgain_outputs15_Jun_2020'; % copo, kmin=0.00001
filename = 'estim_LOMgain_outputs15_Jun_2020_15_45_21'; % copo, kmin=0
load([filename,'.mat'])
xxgrid = estim_outputs{1};
yygrid = estim_outputs{2};
ng     = estim_outputs{3};
k1_opt = estim_outputs{4};
alph_opt = estim_outputs{5};
x = estim_outputs{6};
tol    = estim_outputs{7};
lbname = estim_outputs{8};
ndrop_est  = estim_outputs{9};

alph = alph_opt;

% zoom in on the grid
k1 = ndim_simplex_eval(x,[xxgrid(:)';yygrid(:)'./10],alph_opt);
figure
set(gcf,'color','w'); % sets white background color
set(gcf, 'Position', get(0, 'Screensize')); % sets the figure fullscreen
surf(xxgrid, yygrid./10,reshape(k1,[ng,ng]))
xlabel('$k^{-1}_{t-1}$','interpreter', 'latex', 'fontsize', fs)
ylabel('$fe_{t|t-1}$','interpreter', 'latex', 'fontsize', fs)
zlabel('$k^{-1}_{t}$','interpreter', 'latex', 'fontsize', fs, 'rotation',0)
% title('Optimal','interpreter', 'latex', 'fontsize', fs)
ax = gca; % current axes
ax.FontSize = fs;
set(gca,'TickLabelInterpreter', 'latex');
grid on
grid minor
close

% Specify info assumption on the Taylor rule and not to include a monpol
% shock
knowTR =0
mpshock=0
%%%%%%%%%%%%%%%%%%%

[param, setp, param_names, param_values_str, param_titles] = parameters_next;
[x0, y0, k0, phi0, FA0, FB0, FEt_10, diff0] = sim_learnLH_clean(param,gx,hx,eta,PLM, gain, T,ndrop,e,knowTR,mpshock);
create_plot_observables(y0,seriesnames, 'Simulation using the Taylor rule', [this_code, '_implement_anchTC_obs_TR_',PLM_name,'_', todays_date], print_figs)
create_plot_observables(1./k0,invgain, 'Simulation using the Taylor rule', [this_code, '_implement_anchTC_invgain_TR_',PLM_name,'_', todays_date], print_figs)


return
%% Let's try to fit an AR(1) LOM gain to this estimated one
ARLOM = @(param,k1t_1,fe) param(1).*k1t_1 +param(2).*fe.^2;

param0 = [0.1,0.0001];
k1AR0 = ARLOM(param0,xxgrid,yygrid);

fitAR = @(param,k1t_1,fe, k1_target) (k1_target - vec(ARLOM(param,k1t_1,fe)))'*(k1_target - vec(ARLOM(param,k1t_1,fe)));
loss0 = fitAR(param0,xxgrid,yygrid,k1_opt);

%Optimization Parameters
options = optimset('fmincon');
options = optimset(options, 'TolFun', 1e-9, 'display', 'iter'); 

ub = [1,0.5];
lb = [0,0];
tic
%Declare a function handle for optimization problem
objh = @(param) fitAR(param,xxgrid,yygrid,k1_opt);
[param_opt, loss_opt] = fmincon(objh, param0+1*rand(size(param0)), [],[],[],[],lb,ub,[],options);
% fmincon(FUN,X0,A,B,Aeq,Beq,LB,UB,NONLCON,OPTIONS)
toc

param_opt

k1AR_opt = ARLOM(param_opt,xxgrid,yygrid);
figure
set(gcf,'color','w'); % sets white background color
set(gcf, 'Position', get(0, 'Screensize')); % sets the figure fullscreen
surf(xxgrid, yygrid,reshape(k1AR_opt,[ng,ng]))
xlabel('$k^{-1}_{t-1}$','interpreter', 'latex', 'fontsize', fs)
ylabel('$fe_{t|t-1}$','interpreter', 'latex', 'fontsize', fs)
zlabel('$k^{-1}_{t}$','interpreter', 'latex', 'fontsize', fs, 'rotation',0)
% title('Optimal','interpreter', 'latex', 'fontsize', fs)
ax = gca; % current axes
ax.FontSize = fs;
set(gca,'TickLabelInterpreter', 'latex');
grid on
grid minor
figname = [this_code, '_fitAR_', todays_date];
if print_figs ==1
    disp(figname)
    cd(figpath)
    export_fig(figname)
    cd(current_dir)
    close
end