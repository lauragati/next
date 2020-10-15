% command_IRFs_approx_pretty_together.m
% adapted from command_IRFs_approx_pretty to put all IRFs on the same graph
% 15 Oct 2020

clearvars
close all
clc

date_today = strrep(datestr(today, 'yyyy-mm-dd'), '-','_');

% % Add all the relevant paths and grab the codename
this_code = mfilename;
[current_dir, basepath, BC_researchpath,toolpath,export_figpath,figpath,tablepath,datapath] = add_paths;

% Variable stuff ---
print_figs        = 1;
stop_before_plots = 0;
skip_old_plots    = 0;
output_table = print_figs;

plot_IRFs=0;
plot_simulated_sequence = 0;
plot_gain_IRF = 0;
plot_IRFs_anch = 1; % conditional on being anchored when shock hits, not trivial for smooth anchoring
plot_gains=1;

%% Parameters
tic
[param, setp, param_names, param_values_str, param_titles] = parameters_next;


psi_x = param.psi_x;
psi_pi = param.psi_pi;
lamx =param.lamx;
lami = param.lami;

sig_r = param.sig_r;
sig_i = param.sig_i;
sig_u = param.sig_u;
ne = 3;

% Params for the general learning code
constant_only = 1; % learning constant only
constant_only_pi_only = 11; % learning constant only, inflation only
mean_only_PLM = -1;
slope_and_constant = 2;

dgain = 1;  % 1 = decreasing gain, 2 = endogenous gain, 3 = constant gain
again_critCEMP  = 21;
again_critCUSUM = 22;
again_critsmooth = 23;
cgain = 3;

% filename = 'estim_LOMgain_outputs_univariate_coax11_Sep_2020_15_46_40'; % materials44 candidate
filename = 'estim_LOMgain_outputs_univariate_coax15_Sep_2020_16_14_00'; % complete materials44 candidate (21 Sept draft)

% load the saved stuff
load([filename,'.mat'])
% Structure of saved file:
% estim_configs={nfe,gridspacing,femax,femin,ub,lb,Wprior,Wdiffs2,Wmid,Wmean,T,ndrop,N,eN, rngsetting};
% learn_configs = {param,PLM_name, gain_name, knowTR, mpshock};
% estim_outputs = {fegrid_fine, ng_fine, alph_opt, alph_k, ALPH0, x, estim_configs, learn_configs};
fegrid_fine = estim_outputs{1};
ng_fine     = estim_outputs{2};
alph_opt      = estim_outputs{3};
alph_k        = estim_outputs{4};
ALPH0         = estim_outputs{5};
x             = estim_outputs{6};
estim_configs = estim_outputs{7};
learn_configs = estim_outputs{8};
nfe            = estim_configs{1};
gridspacing    = estim_configs{2};
femax          = estim_configs{3};
femin          = estim_configs{4};
ub             = estim_configs{5};
lb             = estim_configs{6};
Wprior         = estim_configs{7};
Wdiffs2        = estim_configs{8};
Wmid           = estim_configs{9};
Wmean          = estim_configs{10};
T_est          = estim_configs{11};
ndrop_est      = estim_configs{12};
N_est          = estim_configs{13};
eN_est         = estim_configs{14};
rngsetting_est = estim_configs{15};
param       = learn_configs{1};
PLM_name    = learn_configs{2};
gain_name   = learn_configs{3};
knowTR_est  = learn_configs{4};
mpshock_est = learn_configs{5};

% return
alph=alph_opt;
fegrid = x{1};

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% % this here is the 'truth' in acf_sim_univariate_data_06_Jul_2020.
% nfe=5;
% femax = 2;
% femin = -femax;
% fegrid = linspace(femin,femax,nfe);
% x = cell(1,1);
% x{1} = fegrid;
% alph_true = [0.05;0.025;0;0.025;0.05];
% alph=alph_true;

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% %% 27 August 2020: calibration C (Materials 43)
%
% alph = [0.8    0.4         0    0.4    0.8]'
% fegrid = [-4,-3,0,3,4]
% x{1} = fegrid;
%
% [param, setp, param_names, param_values_str, param_titles] = parameters_next;

%% Varying stuff

% 2 std-dev monpol shock
sig_i = 2*sig_i;
param.sig_i = sig_i;

psi_pi_vals = [1.01, 1.5, 2];
P = length(psi_pi_vals);


% Model selection
%%%%%%%%%%%%%%%%%%%
PLM = constant_only_pi_only;
gain = again_critsmooth;
%%%%%%%%%%%%%%%%%%%

T = 400 % 400
% Size of cross-section
N = 100 %100, 500, 1000

knowTR=1
mpshock=1

ndrop = 5 % 5, 100, 0
dt_vals = 25; %25 time of imposing innovation
h = 10; % h-period IRFs

% preallocate IRFs for all values of psi_pi
iry_P = zeros(ne,h,P);
RIR_anch_P = zeros(ne,h,P);
RIR_unanch_P = zeros(ne,h,P);


for p=1:P
    param.psi_pi = psi_pi_vals(p); %1.01, 2
    
    
    
    % RE model
    [fyn, fxn, fypn, fxpn] = model_NK(param);
    [gx,hx]=gx_hx_alt(fyn,fxn,fypn,fxpn);
    [ny, nx] = size(gx);
    [Aa, Ab, As] = matrices_A_13_true_baseline(param, hx);
    SIG = eye(nx).*[sig_r, sig_i, sig_u]';
    eta = SIG; %just so you know
    
    [PLM_name, gain_name, gain_title] = give_names(PLM, gain);
    
    
    %% Simulate models
    
    % gen all the N sequences of shocks at once.
    rng(0)
    eN = randn(ne,T+ndrop,N);
    v = zeros(ne+1, T+ndrop); % shut off measurement error
    
    % Preallocate
    nd = size(dt_vals,2);
    d = 1; % the innovation, delta
    GIR_Y_LH = zeros(ny,h,N,nd);
    GIR_k = zeros(h,N,nd);
    k = zeros(T,N);
    ks= zeros(T,N);
    k_dt = zeros(1,N); % gain when the shock hit
    
    % warning off
    for s=2  %2->zoom in on monetary policy shock
        x0 = zeros(1,nx);
        x0(s) = d;
        
        for n=1:N
            % Sequence of innovations
            e = squeeze(eN(:,:,n));
            
            % Unshocked
            % RE
            [x_RE, y_RE] = sim_model(gx,hx,SIG,T,ndrop,e);
            % Learning
            %         [x_LH, y_LH, ~, ~, ~, ~, ~, ~, diff,~, k(:,n)] = sim_learnLH(gx,hx,SIG,T+ndrop,ndrop,e, Aa, Ab, As, param, PLM, gain);
            if gain ~= again_critsmooth
                [x_LH, y_LH, k(:,n), ~, ~, ~, ~, diff] = sim_learnLH_clean(param,gx,hx,eta, PLM, gain, T+ndrop,ndrop,e,knowTR,mpshock);
            else
                [x_LH, y_LH, k(:,n), ~, ~, ~, ~,diff] = sim_learnLH_clean_approx_univariate(alph,x,param,gx,hx,eta, PLM, gain,T+ndrop,ndrop,e,v,knowTR,mpshock);
            end
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
                %             [~, ys_LH, ~, ~, ~, ~, ~, ~, ~,~, ks(:,n), anch(n)] = sim_learnLH(gx,hx,SIG,T+ndrop,ndrop,e, Aa, Ab, As, param, PLM, gain, dt, x0);
                if gain ~= again_critsmooth
                    [~, ys_LH, ks(:,n), ~, ~, ~, ~, diffs] = sim_learnLH_clean(param,gx,hx,eta, PLM, gain, T+ndrop,ndrop,e,knowTR,mpshock, dt, x0);
                else
                    [~, ys_LH, ks(:,n), ~, ~, ~, ~,diffs, ~,~,~,~,k_dt(n)] = sim_learnLH_clean_approx_univariate(alph,x,param,gx,hx,eta, PLM, gain,T+ndrop,ndrop,e,v,knowTR,mpshock,dt,x0);
                end
                % Construct GIRs
                GIR_Y_LH(:,:,n,t) = ys_LH(:,dt:dt+h-1) - y_LH(:,dt:dt+h-1);
                GIR_k(:,n,t) = ks(dt:dt+h-1,n) - k(dt:dt+h-1,n);
            end
            
        end
    end
    
    % Annualize inflation and interest rates
    GIR_Y_LH([1,3],:,:) =  ((GIR_Y_LH([1,3],:,:)/100+1).^4 -1)*100;
    iry([1,3],:,:) =  ((iry([1,3],:,:)/100+1).^4 -1)*100;
    
    
    % Gather the gains when the shock hit and calculate 10 and 90 percentile in
    % the cross-section
    k1_dt_sort = sort(1./k_dt);
    if mod(N,10)~=0
        k1_10 = k1_dt_sort(floor(N/10)+1);
        k1_90 = k1_dt_sort(ceil(9*N/10)+1);
    elseif mod(N,10)==0
        k1_10 = (k1_dt_sort(floor(N/10)) + k1_dt_sort(floor(N/10)+1)) /2;
        k1_90 = (k1_dt_sort(ceil(9*N/10)) + k1_dt_sort(ceil(9*N/10)+1)) /2;
    end
    well_anch_idx = find(1./k_dt <= k1_10);
    unanch_idx = find(1./k_dt >= k1_90);
    
    % warning on
    % Construct RIRs by simple method: means (Option 1)
    RIR_Y_LH = squeeze(mean(GIR_Y_LH,3));
    RIR_anch = squeeze(mean(GIR_Y_LH(:,:,well_anch_idx),3));
    RIR_unanch = squeeze(mean(GIR_Y_LH(:,:,unanch_idx),3));
    RIR_k = squeeze(mean(GIR_k,2));
    RIR_kinv = RIR_k;
    % only invert for nonzero elements
    RIR_kinv(abs(RIR_k) > 0) = 1./RIR_k(abs(RIR_k) > 0);
    
    disp(['(psi_x, psi_pi, lamx, lami)=   ', num2str([psi_x, psi_pi, lamx, lami])])
    toc
    
    iry_P(:,:,p) = iry;
    RIR_anch_P(:,:,p) = RIR_anch;
    RIR_unanch_P(:,:,p) = RIR_unanch;
    
    
    if stop_before_plots==1
        return
    end
end
% return
%% Plots

shocknames = {'natrate', 'monpol','costpush'};
% titles_obs = {'Inflation','Output gap','Interest rate'};
titles_obs = {'$\pi$','$x$','$i$'};

figspecs = ['_',this_code, '_', date_today];

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%  IRF: OBSERVABLES LH against RE, anchored AND unanchored
series1 = iry_P(:,:,1);
series2 = iry_P(:,:,2);
series3 = iry_P(:,:,3);
series4 = RIR_unanch_P(:,:,1);
series5 = RIR_unanch_P(:,:,2);
series6 = RIR_unanch_P(:,:,3);
subplot_names = titles_obs;
legendnames = {'$\psi_{\pi} = 1.01$', '$\psi_{\pi} = 1.5$', '$\psi_{\pi} = 2$'};
figtitle = [gain_title, '; when shock imposed at t=', num2str(dt_vals(t)), ', unanchored'];

% Plot configs
[fs, lw] = plot_configs;

figure
set(gcf,'color','w'); % sets white background color
set(gcf, 'Position', get(0, 'Screensize')); % sets the figure fullscreen
for i=1:ny
    subplot(1,ny,i)
    plot(zeros(1,length(series1)), 'k', 'linewidth',lw); hold on
%     h1(i) =   plot(series1(i,:), 'color', 'r','linewidth', lw);
%     h2(i) =   plot(series2(i,:), 'color', 'r', 'linewidth', lw, 'linestyle', '--');
%     h3(i) =   plot(series3(i,:), 'color', 'r', 'linewidth', lw, 'marker', 'o');
    h4(i) =   plot(series4(i,:), 'linewidth', lw);
    h5(i) =   plot(series5(i,:), 'linewidth', lw, 'linestyle', '--');
    h6(i) =   plot(series6(i,:), 'linewidth', lw, 'marker', 'o');
    ax = gca; % current axes
    ax.FontSize = fs/1.2;
    set(gca,'TickLabelInterpreter', 'latex');
    grid on
    grid minor
    title(titles_obs{i}, 'interpreter', 'latex', 'fontsize', fs)
    ax.YAxis.Exponent = 0;
    ax.YRuler.Exponent = 0; % turns off scientific notation
end

% Overall legend
% add a bit space to the figure
fig = gcf;
fig.Position(3) = fig.Position(3) + 250;
% add legend
Lgnd = legend('show',[h4(i), h5(i), h6(i)], legendnames, 'location', 'southoutside', 'interpreter', 'latex', 'NumColumns', 3);
legend('boxoff')
Lgnd.Position(1) = 0.38;
Lgnd.Position(2) = 0;

figname = ['RIR_together_psi_pi' ,strrep(num2str(param.psi_pi), '.','_'),figspecs];
if print_figs ==1
    disp(figname)
    cd(figpath)
    export_fig(figname)
    cd(current_dir)
    close
end
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
