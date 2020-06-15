%clean_up_sim_learn
% Clean up sim_learnLH.m - this here copies % command_IRFs_anchoring.m
% The nice thing is you can run this and perform comparisons with the
% retired sim_learnLH.m, and its equivalent cleaned-up version,
% sim_learnLH_clean.m. The last version, sim_learnLH_clean_smooth is for
% future use: it implements the smooth anchoring function version of the
% model. 
% 10 April 2020

clearvars
close all
clc

date_today = strrep(datestr(today, 'yyyy-mm-dd'), '-','_');

% % Add all the relevant paths and grab the codename
this_code = mfilename;
[current_dir, basepath, BC_researchpath,toolpath,export_figpath,figpath,tablepath,datapath,tryouts_path,inputsRyan_path] = add_paths;

% Variable stuff ---
print_figs        = 0;
stop_before_plots = 0;
skip_old_plots    = 0;
output_table = print_figs;

plot_IRFs=1;
plot_simulated_sequence = 1;
plot_gains=0;
plot_gain_IRF = 0;
plot_IRFs_anch = 0;
skip_old_stuff = 1;

% Parameters
tic
[param, set, param_names, param_values_str, param_titles] = parameters_next;


psi_x = param.psi_x;
psi_pi = param.psi_pi;
thetbar =param.thetbar;
thettilde = param.thettilde;

sig_r = param.sig_r;
sig_i = param.sig_i;
sig_u = param.sig_u;
ne = 3;

% Params for the general learning code
constant_only = 1; % learning constant only (vector - _smooth won't work in this case)
constant_only_pi_only = 11; % learning constant only, inflation only (scalar)
slope_and_constant = 2;

dgain = 1;  % 1 = decreasing gain, 2 = endogenous gain, 3 = constant gain
again_critCEMP  = 21;
again_critCUSUM = 22;
again_critsmooth = 23;
cgain = 3;


% Model selection and informational assumption

% Model selection
%%%%%%%%%%%%%%%%%%%
PLM = constant_only_pi_only;
gain = again_critsmooth;
%%%%%%%%%%%%%%%%%%%
% % the smooth criterion only works in this combo:
% PLM = constant_only_pi_only;
% gain = again_critsmooth;

T = 400 % 400
% Size of cross-section
N = 1 %100 500
burnin = 100; %100
dt_vals = 25; %25 time of imposing innovation 345
h = 10; % h-period IRFs

% RE model
[fyn, fxn, fypn, fxpn] = model_NK(param);
[gx,hx]=gx_hx_alt(fyn,fxn,fypn,fxpn);
[ny, nx] = size(gx);
[Aa_old, Ab_old, As_old] = matrices_A_13_true_baseline(param, hx);
[Aa_april, Ab_april, As_april] = matrices_A_25_true_baseline(param,hx);

Aa = Aa_april;
Ab = Ab_april;
As = As_april;

SIG = eye(nx).*[sig_r, sig_i, sig_u]';
eta = SIG; %just so you know

[PLM_name, gain_name, gain_title] = give_names(PLM, gain);


% Simulate models

% gen all the N sequences of shocks at once.
rng(0)
eN = randn(ne,T+burnin,N);

% Preallocate
nd = size(dt_vals,2);
d = 1; % the innovation, delta
GIR_Y_LH = zeros(ny,h,N,nd);
GIR_k = zeros(h,N,nd);
k = zeros(T,N);
ks= zeros(T,N);
k_clean = zeros(T,N);
k_smooth = zeros(T,N);
anch = zeros(1,N); % indicator for whether that sequence was anchored when the shock hit
g_pi = zeros(T,1);
g_pibar = zeros(T,1);
% warning off
for s=2  %2->zoom in on monetary policy shock
    x0 = zeros(1,nx);
    x0(s) = d;
    
    for n=1:N
        % Sequence of innovations
        e = squeeze(eN(:,:,n));
        
        % Unshocked
        %         dbstop if warning
        % RE
        [x_RE, y_RE] = sim_model(gx,hx,SIG,T,burnin,e);
        % Learning
        % the original code
        [x_LH, y_LH, ~, ~, ~, ~, ~, ~, diff(:,n),~, k(:,n)] = sim_learnLH(gx,hx,SIG,T+burnin,burnin,e, Aa, Ab, As, param, PLM, gain);
        % its plain vanilla cleaned version WORK WITH THIS
        [x_clean, y_clean, k_clean(:,n), phi_clean, ~, ~, ~,diff_clean(:,n)] = sim_learnLH_clean(param,gx,hx,SIG,PLM, gain, T+burnin,burnin,e);
        
        % RETIRED
%         % a cleaned version only for pi-only scalar learning with scalar
%         % smooth criterion RETIRED
%         [x_smooth, y_smooth, k_smooth(:,n), pibar(:,n), FA, FB, g_pi, g_pibar, ~,diff_smooth(:,n)] = sim_learnLH_clean_smooth(param,gx,hx,SIG, Aa, Ab, As,T+burnin,burnin,e);
%         % vector learning for the smooth anchoring function g. RETIRED
%         [x_g, y_g, k_g(:,n),  phi_g, ~, ~, diff_g(:,n)] = sim_learnLH_clean_g(param,gx,hx,SIG, Aa, Ab, As,PLM, T+burnin,burnin,e);
        % Shocked
        % RE
        % make RE shock the same scale as learning:
        x0RE = (SIG*x0')';
        [IR, iry, irx]=ir(gx,hx,x0RE,h);
        iry = iry';
        RE_fcsts = gx*hx*irx';
        
        % Learning
        for t=1:nd
            dt = dt_vals(t);
            
            % Shocked
            [~, ys_LH] = sim_learnLH(gx,hx,SIG,T+burnin,burnin,e, Aa, Ab, As, param, PLM, gain, dt, x0);
            [~, ys_clean] = sim_learnLH_clean(param,gx,hx,SIG,PLM, gain, T+burnin,burnin,e, dt, x0);
            % RETIRED:
%             [~, ys_smooth] = sim_learnLH_clean_smooth(param,gx,hx,SIG, Aa, Ab, As,T+burnin,burnin,e,dt,x0);
%             [~, ys_g] = sim_learnLH_clean_g(param,gx,hx,SIG, Aa, Ab, As,PLM, T+burnin,burnin,e, dt, x0);
        end
        
    end
end

disp('%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%')
disp('%%%%%  Diffs for sim_learnLH_clean.m %%%%%')
disp('%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%')
disp(['Max difference in states = ', num2str(max(max(abs(x_LH-x_clean))))])
disp(['Max difference in jumps = ', num2str(max(max(abs(y_LH-y_clean))))])
disp(['Max difference in k = ', num2str(max(max(abs(k-k_clean))))])
disp(['Max difference in convergence diffs = ', num2str(max(max(abs(diff-diff_clean))))])
disp(['Max difference in shocked jumps = ', num2str(max(max(abs(ys_LH-ys_clean))))])
% disp('%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%')
% disp('%%%%%  Diffs for sim_learnLH_clean_smooth.m %%%%%')
% disp('%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%')
% disp(['Max difference in states = ', num2str(max(max(abs(x_LH-x_smooth))))])
% disp(['Max difference in jumps = ', num2str(max(max(abs(y_LH-y_smooth))))])
% disp(['Max difference in k = ', num2str(max(max(abs(k-k_smooth))))])
% disp(['Max difference in convergence diffs = ', num2str(max(max(abs(diff-diff_smooth))))])
% disp(['Max difference in shocked jumps = ', num2str(max(max(abs(ys_LH-ys_smooth))))])
% disp('%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%')
% disp('%%%%% Diffs for sim_learnLH_clean_g.m %%%%%')
% disp('%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%')
% disp(['Max difference in states = ', num2str(max(max(abs(x_LH-x_g))))])
% disp(['Max difference in jumps = ', num2str(max(max(abs(y_LH-y_g))))])
% disp(['Max difference in k = ', num2str(max(max(abs(k-k_g))))])
% disp(['Max difference in convergence diffs = ', num2str(max(max(abs(diff-diff_g))))])
% disp(['Max difference in shocked jumps = ', num2str(max(max(abs(ys_LH-ys_g))))])


%%
figure
subplot(1,2,1)
plot(mean(diff_clean,2))
title('Convergence - mean(diff) across N')
subplot(1,2,2)
plot(diff_clean(:,end))
title('Convergence  - diff for N=end')
