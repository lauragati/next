% command_check_simulation_approx.m
% a code meant to see how simulating the model reacts to various alphas.
% also look at how various simulations affect the estimation loss Nsimul.
% 10 August 2020

clearvars
close all
clc

% Add all the relevant paths and grab the codename
this_code = mfilename;
[current_dir, basepath, BC_researchpath,toolpath,export_figpath,figpath,tablepath,datapath, inputsRyan_path] = add_paths;
todays_date = strrep(datestr(today), '-','_');
nowstr = strrep(strrep(strrep(datestr(now), '-','_'), ' ', '_'), ':', '_');

% Variable stuff ---
print_figs        = 0;
stop_before_plots = 0;
skip_old_plots    = 0;
output_table = print_figs;

skip = 0;
[fs, lw] = plot_configs;
datestr(now)


save_stuff=0;


%% Params, initialization

% Grid
nfe = 5 % 5,7,9
gridspacing = 'uniform'; % uniform or uneven
% grids for fe_{t|t-1}
femax = 2; % 3.5
femin = -femax;
sig_v=0;
nobs=3;

T=155
N=100
ndrop = 5

rng(1) % rng(1)  vs. rng('default')=rng(0) is the one that was used to generate the true data.

% gen all the N sequences of shocks at once.
eN = randn(3,T+ndrop,N);
vN = sig_v*randn(nobs,T+ndrop,N); % measurement error


[param, setp, param_names, param_values_str, param_titles] = parameters_next;

sig_r = param.sig_r;
sig_i = param.sig_i;
sig_u = param.sig_u;

% RE model
[fyn, fxn, fypn, fxpn] = model_NK(param);
[gx,hx]=gx_hx_alt(fyn,fxn,fypn,fxpn);
[~, nx] = size(gx);

% [Aa, Ab, As] = matrices_A_13_true_baseline(param, hx);
SIG = eye(nx).*[sig_r, sig_i, sig_u]';
eta = SIG; %just so you know

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

% Model selection
%%%%%%%%%%%%%%%%%%%
PLM = constant_only_pi_only;
gain = again_critsmooth;
%%%%%%%%%%%%%%%%%%%

% display which model you're doing
[PLM_name, gain_name, gain_title] = give_names(PLM, gain);

% Specify info assumption on the Taylor rule and whether to include a monpol
% shock
knowTR =1
mpshock=1
%%%%%%%%%%%%%%%%%%%
% turned monpol shocks on in smat.m to avoid stochastic singularity!

switch gridspacing
    case 'uniform'
        fegrid = linspace(femin,femax,nfe); % for alph0, fe is between (-2.6278,3.5811).
    case 'uneven'
        fegrid = uneven_grid(femin,femax,nfe)
    case 'manual'
        fegrid = [-2,-1,1,2]
end
% map to ndim_simplex
x = cell(1,1);
x{1} = fegrid;
[xxgrid] = meshgrid(fegrid);


% Finer sample
ng_fine = 100;
broaden = 2
fegrid_fine = linspace(femin-broaden,femax+broaden,ng_fine);

figspecs = ['_',PLM_name,'_', 'N_', num2str(N),'_nfe_', num2str(nfe), '_femax_', num2str(femax), '_', this_code, '_', todays_date];

if skip==0
    
    %% I) Check simulation at truth for various error sequences.
    
    alph_true = 2.5*[0.05;0.025;0;0.025;0.05];
    y = zeros(nobs,T,N);
    k = zeros(T,N);
    pibar = zeros(T,N);
    fe = zeros(T,N);
    diffs = nan(T+ndrop,N);
    
    tic
    parfor n=1:N
        e_n = squeeze(eN(:,:,n));
        v_n = squeeze(vN(:,:,n));
        [~, y(:,:,n), k(:,n),phi,~,~,FE, diffs(:,n)] = sim_learnLH_clean_approx_univariate(alph_true,x,param,gx,hx,eta, PLM, gain, T+ndrop,ndrop,e_n,v_n, knowTR,mpshock);
        pibar(:,n) = squeeze(phi(1,1,:));
        fe(:,n) = FE(1,:);
    end
    toc
    % Cut off last, unset period
    k1 = 1./k(1:end-1,:);
    y = y(:,1:end-1,:);
    pibar = pibar(1:end-1,:);
    fe = fe(1:end-1,:);
    
    
    %% Plots I
    
    figure
    set(gcf,'color','w'); % sets white background color
    set(gcf, 'Position', get(0, 'Screensize')); % sets the figure fullscreen
    subplot(1,3,1)
    histogram(k1)
    ax = gca; % current axes
    ax.FontSize = fs*3/4;
    set(gca,'TickLabelInterpreter', 'latex');
    grid on
    grid minor
    title('$k^{-1}$ in cross-section','interpreter', 'latex', 'fontsize', fs*3/4)
    subplot(1,3,2)
    histogram(pibar)
    ax = gca; % current axes
    ax.FontSize = fs*3/4;
    set(gca,'TickLabelInterpreter', 'latex');
    grid on
    grid minor
    title('$\bar{\pi}$ in cross-section','interpreter', 'latex', 'fontsize', fs*3/4)
    subplot(1,3,3)
    histogram(fe)
    ax = gca; % current axes
    ax.FontSize = fs*3/4;
    set(gca,'TickLabelInterpreter', 'latex');
    grid on
    grid minor
    title('fe in cross-section','interpreter', 'latex', 'fontsize', fs*3/4)
    
    figname = ['distribs', figspecs];
    if print_figs ==1
        disp(figname)
        cd(figpath)
        export_fig(figname)
        cd(current_dir)
        close
    end
    
    figure
    set(gcf,'color','w'); % sets white background color
    set(gcf, 'Position', get(0, 'Screensize')); % sets the figure fullscreen
    subplot(1,3,1)
    plot(mean(k1,2), 'linewidth', lw)
    ax = gca; % current axes
    ax.FontSize = fs*3/4;
    set(gca,'TickLabelInterpreter', 'latex');
    grid on
    grid minor
    title('Cross-sectional mean of $k^{-1}$','interpreter', 'latex', 'fontsize', fs*3/4)
    subplot(1,3,2)
    plot(mean(pibar,2), 'linewidth', lw)
    ax = gca; % current axes
    ax.FontSize = fs*3/4;
    set(gca,'TickLabelInterpreter', 'latex');
    grid on
    grid minor
    ax.YRuler.Exponent = 0; % turns off scientific notation
    title('Cross-sectional mean of $\bar{\pi}$','interpreter', 'latex', 'fontsize', fs*3/4)
    subplot(1,3,3)
    plot(mean(fe,2), 'linewidth', lw)
    ax = gca; % current axes
    ax.FontSize = fs*3/4;
    set(gca,'TickLabelInterpreter', 'latex');
    grid on
    grid minor
    title('Cross-sectional mean of fe','interpreter', 'latex', 'fontsize', fs*3/4)
    
    figname = ['means', figspecs];
    if print_figs ==1
        disp(figname)
        cd(figpath)
        export_fig(figname)
        cd(current_dir)
        close
    end
    
end

return
%% Check cross-sections for various alphas
alph_true = [0.05;0.025;0;0.025;0.05];

scenarios = {[1,5], [2,4], 3};
nscen = numel(scenarios);
nrange =3;
alphi_values =linspace(0,0.1,nrange);

mean_k1 = nan(T-1,nrange,nscen);
mean_pibar = nan(T-1,nrange,nscen);
mean_fe = nan(T-1,nrange,nscen);

tic
for s =1:nscen
    idx = scenarios{s};
    for i=1:nrange
        alph = alph_true;
        alph(idx) = alphi_values(i);
        
        k1 = zeros(T-1,N);
        pibar = zeros(T-1,N);
        fe = zeros(T-1,N);
        parfor n=1:N
            e_n = squeeze(eN(:,:,n));
            v_n = squeeze(vN(:,:,n));
            try
                [~, ~, k, phi,~,~,FE] = sim_learnLH_clean_approx_univariate(alph,x,param,gx,hx,eta, PLM, gain, T+ndrop,ndrop,e_n,v_n, knowTR,mpshock);
                k1(:,n) = 1./k(1:end-1);
                pibar(:,n) = squeeze(phi(1,1,1:end-1));
                fe(:,n) = FE(1,1:end-1);
            catch
            end
        end
        mean_k1(:,i,s) = mean(k1,2);
        mean_pibar(:,i,s) = mean(pibar,2);
        mean_fe(:,i,s) = mean(fe,2);
    end
end
toc

%% For each scenario and each alpha, plot means

% k1
for s=1:nscen
    idx = scenarios{s};
    figure
    set(gcf,'color','w'); % sets white background color
    set(gcf, 'Position', get(0, 'Screensize')); % sets the figure fullscreen
    for i=1:nrange
        alph = alph_true;
        alph(idx) = alphi_values(i);
        
        subplot(nrange,1,i)
        plot(squeeze(mean_k1(:,i,s)), 'linewidth', lw)
        ax = gca; % current axes
        ax.FontSize = fs*3/4;
        set(gca,'TickLabelInterpreter', 'latex');
        grid on
        grid minor
        title(['$\alpha =$  ', num2str(alph')],'interpreter', 'latex', 'fontsize', fs*3/4)
    end
    sgt = sgtitle(['$k^{-1}$ - varying $\alpha$ ', num2str(idx)]);
    sgt.FontSize =fs;
    sgt.Interpreter = 'latex';
    
    figname = ['k1_scenario_', num2str(s), '_', figspecs];
    if print_figs ==1
        disp(figname)
        cd(figpath)
        export_fig(figname)
        cd(current_dir)
        close
    end
end


% pibar
for s=1:nscen
    idx = scenarios{s};
    figure
    set(gcf,'color','w'); % sets white background color
    set(gcf, 'Position', get(0, 'Screensize')); % sets the figure fullscreen
    for i=1:nrange
        alph = alph_true;
        alph(idx) = alphi_values(i);
        
        subplot(nrange,1,i)
        plot(squeeze(mean_pibar(:,i,s)), 'linewidth', lw)
        ax = gca; % current axes
        ax.FontSize = fs*3/4;
        set(gca,'TickLabelInterpreter', 'latex');
        grid on
        grid minor
        title(['$\alpha =$  ', num2str(alph')],'interpreter', 'latex', 'fontsize', fs*3/4)
    end
    sgt = sgtitle(['$\bar{\pi}$ - varying $\alpha$ ', num2str(idx)]);
    sgt.FontSize =fs;
    sgt.Interpreter = 'latex';
    
    figname = ['pibar_scenario_', num2str(s), '_', figspecs];
    if print_figs ==1
        disp(figname)
        cd(figpath)
        export_fig(figname)
        cd(current_dir)
        close
    end
end

% fe
for s=1:nscen
    idx = scenarios{s};
    figure
    set(gcf,'color','w'); % sets white background color
    set(gcf, 'Position', get(0, 'Screensize')); % sets the figure fullscreen
    for i=1:nrange
        alph = alph_true;
        alph(idx) = alphi_values(i);
        
        subplot(nrange,1,i)
        plot(squeeze(mean_fe(:,i,s)), 'linewidth', lw)
        ax = gca; % current axes
        ax.FontSize = fs*3/4;
        set(gca,'TickLabelInterpreter', 'latex');
        grid on
        grid minor
        title(['$\alpha =$  ', num2str(alph')],'interpreter', 'latex', 'fontsize', fs*3/4)
    end
    sgt = sgtitle(['$fe$ - varying $\alpha$ ', num2str(idx)]);
    sgt.FontSize =fs;
    sgt.Interpreter = 'latex';
    
    figname = ['fe_scenario_', num2str(s), '_', figspecs];
    if print_figs ==1
        disp(figname)
        cd(figpath)
        export_fig(figname)
        cd(current_dir)
        close
    end
end