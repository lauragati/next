%% Coax the solver to get to the right answer: for 1D case only
% older code (June?), adapted 10 Sept 2020 to initialize estimation of real data at
% many points and take the one with lowest loss
% No longer relying on fun_GMM_LOMgain_univariate.m

clearvars
close all
clc

% Add all the relevant paths and grab the codename
this_code = mfilename;
[current_dir, basepath, BC_researchpath,toolpath,export_figpath,figpath,tablepath,datapath] = add_paths;
todays_date = strrep(datestr(today), '-','_');
nowstr = strrep(strrep(strrep(datestr(now), '-','_'), ' ', '_'), ':', '_');

% Variable stuff ---
print_figs        = 0;
if contains(current_dir, 'gsfs0') % sirius server
    print_figs=1;
end
stop_before_plots = 0;
skip_old_plots    = 0;
output_table = print_figs;

save_estim_outputs =0;

skip = 1;
[fs, lw] = plot_configs;
datestr(now)



%% This initialization bit comes straight from command_GMM_LOMgain_univariate.m

filename ='acf_data_21_Jul_2020'; % real data with SPF expectation in it


% Grid
nfe = 5 % 5,7,9
gridspacing = 'manual'; % uniform or uneven, or manual
% grids for fe_{t|t-1}
% femax = 2; % 3.5
% femin = -femax;
% upper and lower bounds for estimated coefficients
ub = ones(nfe,1); %1
lb = zeros(nfe,1); %0
% weights on additional moments
Wprior=0;%0
Wdiffs2= 100000;%100000, seems like 100K is sufficient, or even 10K
Wmid =1000; %1000
Wmean=0;%100, 0
% rng(8)
% alph0 = rand(nfe,1);
% alph0 = 0.1*ones(nfe,1);
alph0 = [0.8,0.5,0,0.5,0.8]'
% alph0 = [0.8,0.4,0,0.4,0.8]'
% alph0 = [0.6,0.3,0,0.3,0.6]'
% alph0 = [0.3,0.1,0,0.1,0.3]'

use_smart_alph0=0;% default
cross_section = 'Nsimulations'; % Nestimations or Nsimulations (default)
est_shocks = 0;
if est_shocks==1
    cross_section = 'Nsimulations';
end
scaleW =0; %0
ryans_rescale = 0; % Per Ryan meeting 12 August 2020.
use_expectations_data=1; %1
sig_v = 0; %0 vs 1 variance of measurement error: set to zero to shut measurement error off (default)

%Optimization Parameters
options = optimoptions('lsqnonlin');
options = optimoptions(options, 'display', 'iter');
% options.TolFun= 1e-9; % objective function tolerance, default 1.0000e-06
% options.OptimalityTolerance = 1e-9; % this is the guy you can access in
% optimoptions, not in optimset. It pertains to first order optimality
% measure. Default 1.0000e-06
options.MaxFunEvals = 700; % 700
% options.MaxIter = 33;
% options.TolX = 1e-9; % step tolerance: default 1.0000e-06
options.UseParallel = 0; % 2/3 of the time
h_sig = 100000;
h_alph= 100;%100000
% options.FiniteDifferenceStepSize = sqrt(eps)*[h_alph;h_alph;h_alph;h_alph;h_alph;]; % default is sqrt(eps); sqrt(eps)*[h_sig;h_sig;h_sig;h_alph;h_alph;h_alph;h_alph;h_alph;]
%%%%%%%%%%%%%%%%%%%

load([filename, '.mat'])
Om_data = acf_outputs{1}; % this is the moments vector
Om_boot = acf_outputs{2}; % moments vectors in bootstrapped samples
nobs = acf_outputs{3};
p = acf_outputs{4};
K = acf_outputs{5};
filt_data = acf_outputs{6};
lost_periods = acf_outputs{7};
T = length(filt_data);
T = T+lost_periods; % to make up for the loss due to filtering

ndrop = 5 % 0-50


% Size of cross-section
N=100

if use_expectations_data==0 % take out moments pertaining to expectations
    Ommatrix = reshape(Om_data,nobs,nobs,K+1);
    % do the same for the bootstrapped variances
    nboot=size(Om_boot,2);
    Om_boot_matrix = reshape(Om_boot, nobs,nobs,K+1, nboot);
    nobs = 3;
    Ommatrix_net = Ommatrix(1:nobs,1:nobs,:);
    Om_boot_matrix_net = Om_boot_matrix(1:nobs,1:nobs,:,:);
    Om_data = vec(Ommatrix_net);
    [omb_m, omb_n, ~, ~] = size(Om_boot_matrix_net);
    Om_boot = reshape(Om_boot_matrix_net,omb_m*omb_n*(K+1), nboot);
    
end

% return

rng(1) % rng(1)  vs. rng('default')=rng(0) is the one that was used to generate the true data.


% gen all the N sequences of shocks at once.
eN = randn(3,T+ndrop,N);
vN = sig_v*randn(nobs,T+ndrop,N); % measurement error

% return

if contains(filename,'sim')
    alph_true = acf_outputs{8};
    nfe_true  = acf_outputs{9};
end

% Note: 3 moments at lag 0 are repeated. So technically we only have
% n_moment-3 moments (and they could be correlated further)

% return

% weighting matrix for GMM
% note: 2nd argument of var(X,W,DIM) signifies whether to normalize by N-1
% or N. 0 -> N-1, 1 -> N. The third argument, DIM, says along which
% dimension to take the variance.
sigboot = diag(var(Om_boot,0,2));
W = sigboot;

% % % Check where the not-rescaling or rescaling causes problems
% X=W;
% scaler = floor(log10(min(diag(W))));
% % a=10^(abs(scaler));
% a=10;
% Y = a*X;
% 1./Y == (1/a) ./X % aha! Not always true!
% inv(Y) == 1/a *inv(X) % but it's b/c 1/a is too close to zero.
%
% 1/a *inv(X) == inv(a*X) % not always true either! So this means that inverting the unrescaled W indeed caused troubles!
% 1/a ./X == 1./(a*X)
% inv(Y) == inv(a*X) % this is always true: so once you scale X, inverting works the way it should. But not before.
%
% inv(X)./inv(Y)
% inv(X)/inv(Y)
%
% % '/' is much more accurate than inv() or ^(-1)!
% X = [1,2; 3,4]
% inv(X) % correct
% 1./X % incorrect wtf!
% inv(X) == 1./X % aren't equal!!! Of course not. Idiot.

% If W really small, scale it up by the exponent of the smallest value
% get exponent of 10 in the smallest diagonal element
% return
if scaleW==1
    scaler = floor(log10(min(diag(sigboot))));
    if scaler < 0
        W_alt = W.* 10^(abs(scaler)); % just to check elementwise, but it gives the same thing
        W = W* 10^(abs(scaler));
    end
end
W1 = W^(-1);
% W1 = W1*10^(abs(scaler)); % scale it back up - do we get what we would w/o scaling? YES.
% W1 - inv(sigboot) % equal on the order of e-12

% % Just checking that it's true for W1 too, and it is
% W1 == 1/a *inv(X) % not always tru
% W1 == inv(a*X) % this is always true
% W1=eye(size(W1));

% Per Ryan meeting 12 August 2020:
W1 = sqrt(W1); % elementwise sqrt.

% Also per Ryan meeting 12 August 2020:
diag_sigboot = diag(sigboot);
maxSig = max(max(diag_sigboot));
minSig = min(min(diag_sigboot));
maxminratio = maxSig/minSig; % Should be and is < 10e+6 if no expectations are used, on the order of 9e+7 or 1e+8 if expectations are used
%  is < 10e+6 if expectations are annualized too in the true data!

if ryans_rescale==1
    % Ryan's rescaling strategy
    threshold_ratio = 10e6;
    sigdiffs = maxSig ./ diag_sigboot;
    % how many too small elements do we have
    n_problems = sum(sigdiffs > threshold_ratio);
    % find the n problematic elements
    where = find(sigdiffs > threshold_ratio, n_problems, 'first');
    % fix the problematic elements of sigma to be 1/10000 the largest element
    % of sigma.
    diag_sigboot(where) = maxSig/10000;
    W = diag(diag_sigboot);
    W1 = inv(W);
end

% return

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
        %         fegrid = [-2,-1,1,2]
        %         fegrid = [-2,-1.5,0,1.5,2]
        fegrid = [-4,-3,0,3,4]
        %         fegrid = [-4,-3, 3,4]
end
% map to ndim_simplex
x = cell(1,1);
x{1} = fegrid;
[xxgrid] = meshgrid(fegrid);

% return
% values for k^{-1}_t for the grid
param.rho_k =0;
kmesh = fk_smooth_pi_only(param,xxgrid,rand(size(xxgrid))); % I've checked and this gives the same as putting fe and k_{t-1} one-by-one thru fk_smooth_pi_only
k1 = 1./kmesh;

% Do an initial approx of the anchoring function to initialize the coeffs
if use_smart_alph0==1
    alph0 = ndim_simplex(x,xxgrid(:)',k1);
end

femin = min(fegrid);
femax = max(fegrid);
% Finer sample
ng_fine = 100;
broaden = 2
fegrid_fine = linspace(femin-broaden,femax+broaden,ng_fine);
k10 = ndim_simplex_eval(x,fegrid_fine,alph0);

% % extrapolation works beautifully
% ndim_simplex_eval(x,[-3:1:3],alph0)


% % Can recover scaling factor?
% res0 = obj_GMM_LOMgain_univariate_mean(alph_true,x,fegrid_fine,param,gx,hx,eta,eN,vN,T,ndrop,PLM,gain,p,Om,W1,Wdiffs2,Wmid,Wmean,use_expectations_data,N);
% resnorm0 = sum(res0.^2);
% % return
% scaler = floor(log10(min(diag(sigboot))));
% sclf = 10^(abs(scaler));
% res1 = obj_GMM_LOMgain_univariate_mean(alph_true,x,fegrid_fine,param,gx,hx,eta,eN,vN,T,ndrop,PLM,gain,p,Om,sqrt(inv(W*sclf)),Wdiffs2,Wmid,Wmean,use_expectations_data,N);
% resnorm1 = sum(res1.^2);

% resnorm0 - resnorm1*sclf
% resnorm0 - resnorm1*sclf^2 % yes you do


% return

%% Pick nsearch starting values

nsearch = 10;


%Declare a function handle for optimization problem
objh = @(alph) obj_GMM_LOMgain_univariate_mean(alph,x,fegrid_fine,param,gx,hx,eta,eN,vN,T,ndrop,PLM,gain,p,Om_data,W1,Wdiffs2,Wmid,Wmean,use_expectations_data,N);
% Try it out once
tic
[res0, Om0, FE0, Om_n0] = obj_GMM_LOMgain_univariate_mean(alph0,x,fegrid_fine,param,gx,hx,eta,eN,vN,T,ndrop,PLM,gain,p,Om_data,W1,Wdiffs2,Wmid,Wmean,use_expectations_data,N);
toc
resnorm0 = sum(res0.^2)


alph1 = ones(nfe,nsearch);
res1 = zeros(length(res0),nsearch);
Om1 = zeros(length(Om_data),nsearch);
flag = zeros(1,nsearch);

% Convex random starting values
if nsearch ==1
    ALPH0 = alph0;
elseif nsearch > 1
    ALPH0 = alph0.^(1:nsearch);
end

start_search = datetime('now');
for i=1:nsearch
    alph0 = ALPH0(:,i);
    %     [alph_opt(:,i), resnorm(i),res(:,i), Om_opt(:,i),flag(i)] = fun_GMM_LOMgain_univariate(acf_outputs, nfe, k1min, k1max, femin, femax, alph0); % retired
    
    tic
    [alph1(:,i), ~,~,flag(i)] = lsqnonlin(objh,alph0,lb,ub,options);
    disp(['Finished ', num2str(i), 'th search'])
    toc
    
    [res1(:,i), Om1(:,i)] = obj_GMM_LOMgain_univariate_mean(alph1(:,i),x,fegrid_fine,param,gx,hx,eta,eN,vN,T,ndrop,PLM,gain,p,Om_data,W1,Wdiffs2,Wmid,Wmean,use_expectations_data,N);
    
end
disp(['Done with all ', num2str(nsearch), ' searches.'])
end_search =datetime('now');
duration_search = end_search - start_search



% Calc sum of squared residuals
loss = sum(res1.^2,1);

disp([num2str(sum(flag>1)/nsearch *100), '% converged'])

% Take the converged ones
alph_conv = alph1(:,flag > 0);
loss_conv = loss(:, flag>0);
Om_conv = Om1(:, flag>0);


% Take the k best points
k=3;
[mink_resnorm,mink_idx]= mink(loss_conv, k);
% and the very best point
[min_resnorm,min_idx]= min(loss_conv);


% The k candidates are
alph_k   = alph_conv(:,mink_idx);
alph_opt = alph_conv(:,min_idx);
Om_k = Om_conv(:,mink_idx);


%% Nice plots of the k best points (also adapted from command_GMM_LOMgain_univariate.m)

% the appendix of each FIGNAME:
figspecs = ['N_', num2str(N),'_nfe_', num2str(nfe), ...
    '_gridspacing_', gridspacing, '_Wdiffs2_', num2str(Wdiffs2),'_Wmid_', num2str(Wmid), '_', cross_section, '_', this_code, '_', nowstr];


% 1.) optimal k alphas
figname = ['alph_opt_',figspecs];
xlab = 'Forecast error';
ylab = 'Gain';
xlplus = [0,0];
ylplus = [0.8,0.4];
legendentries = {};

create_pretty_plot_x_holdon(fegrid,alph_k',legendentries, xlab,ylab,xlplus, ylplus,figname,print_figs)

% return

k = size(Om_k,2);
% 2.) Covariogram
Gamj_data = reshape(Om_data,nobs,nobs,K+1);
Gamj_k = reshape(Om_k,nobs,nobs,K+1,k);

titles = {'$\pi_t$', '$x_t$', '$i_t$', '$E_t\pi_{t+1}$'};
titles_k = {'$\pi_{t-k}$', '$x_{t-k}$', '$i_{t-k}$', '$E_{t-k}\pi_{t-k+1}$'};
it=0;
figure
set(gcf,'color','w'); % sets white background color
set(gcf, 'Position', get(0, 'Screensize')); % sets the figure fullscreen
for i=1:nobs
    for j=1:nobs
        it=it+1;
        sp(it)=subplot(nobs,nobs,it);
        pos_sp = get(sp(it), 'position');
        set(sp(it), 'position', [1, 1, 0.95, 0.95].*pos_sp );
        z = plot(0:K,zeros(1,K+1), 'k--', 'linewidth',lw); hold on
        for kk=1:k
            hk(kk) = plot(0:K,squeeze(Gamj_k(i,j,:,kk)), 'linewidth', lw);
        end
        h = plot(0:K,squeeze(Gamj_data(i,j,:)), 'linewidth', lw);
        
        ax = gca; % current axes
        ax.FontSize = fs*3/4;
        ax.YRuler.Exponent = 0; % turns off scientific notation
        set(gca,'TickLabelInterpreter', 'latex');
        grid on
        grid minor
        title([titles{i}, ' vs. ', titles_k{j}],'interpreter', 'latex', 'fontsize', fs*2/4)
    end
end
% To avoid an undefined property 'NumColumns' on the server:
if contains(current_dir, 'BC_Research') % local
    lh = legend([h, hk(1)],{'Data', 'Best'},'interpreter', 'latex','Position',[0.45 -0.05 0.1 0.2], 'NumColumns',3, 'Box', 'off');
elseif contains(current_dir, 'gsfs0') % sirius server
    lh = legend([h, hk(1)],{'Data', 'Best'},'interpreter', 'latex','Position',[0.45 -0.05 0.1 0.2], 'Box', 'off');
end


% Note position: left, bottom, width, height
figname = ['autocovariogram_' figspecs];
if contains(filename,'sim')==1
    figname = ['autocovariogram_sim_', figspecs];
end
if print_figs ==1
    disp(figname)
    cd(figpath)
    export_fig(figname)
    cd(current_dir)
    close
end

% 3) Histogram of converged losses
figure
histogram(loss_conv)