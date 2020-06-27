% materials35
% Continue to estimate the approximated anchoring function, now focus on
% the 1D case, and initially only on simulated data.
% 24 June 2020

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

skip = 1;
[fs, lw] = plot_configs;
datestr(now)

do21 = 0;
do22 = 0;
do23 = 0;
do24 = 0;
do25 = 1;


%% Hopefully general params

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%% Truth (Simulated data) %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

filename = 'acf_sim_univariate_data_21_Jun_2020'; % simulated data, nfe = 6. Note: I'm using the large moments vector.

load([filename, '.mat'])
Om = acf_outputs{1}; % this is the moments vector
Om_boot = acf_outputs{2}; % moments vectors in bootstrapped samples
ny = acf_outputs{3};
p = acf_outputs{4};
K = acf_outputs{5};
filt_data = acf_outputs{6};
lost_periods = acf_outputs{7};
T = length(filt_data);
T = T+lost_periods; % to make up for the loss due to filtering

if contains(filename,'sim')
    alph_true = acf_outputs{8};
    nfe_true  = acf_outputs{9};
end

% Note: 3 moments at lag 0 are repeated. So technically we only have 42
% moments (and they could be correlated further)
% reshape(Om(1:9),3,3)

% return

% weighting matrix for GMM
% note: 2nd argument of var(X,W,DIM) signifies whether to normalize by N-1
% or N. 0 -> N-1, 1 -> N. The third argument, DIM, says along which
% dimension to take the variance.
W = diag(var(Om_boot,0,2));
W1 = W^(-1);

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%% Model %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
[param, setp, param_names, param_values_str, param_titles] = parameters_next;
ne = 3;

sig_r = param.sig_r;
sig_i = param.sig_i;
sig_u = param.sig_u;

% RE model
[fyn, fxn, fypn, fxpn] = model_NK(param);
[gx,hx]=gx_hx_alt(fyn,fxn,fypn,fxpn);
[ny, nx] = size(gx);

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

% Specify info assumption on the Taylor rule and not to include a monpol
% shock
knowTR =1
mpshock=1
%%%%%%%%%%%%%%%%%%%

% Size of cross-section
% we're not doing a whole cross-section here
ndrop = 5 % 0-50

% gen all the N sequences of shocks at once.
rng(0)
e = randn(ne,T+ndrop);

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%% Approx %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

nfe=6;
k1min = 0;
k1max=1;
femax = 5;
femin = -femax;
fegrid = linspace(femin,femax,nfe);
x = cell(1,1);
x{1} = fegrid;
ng_fine = 100;
fegrid_fine = linspace(femin,femax,ng_fine);

%% 2.1) Evaluate loss on a 6^6 grid
if do21==1
    if skip==0
        mat35_21_server % <--- takes 100 min (1 hr 40 min)
        ngrid = 6;
        alph = zeros(ngrid, 1);
        alphi_values =linspace(0,1,ngrid);
        loss = zeros(ngrid,ngrid);
        tic
        disp('Evaluating 6^6 grid loss, should take 80 min')
        datestr(now)
        for i=1:ngrid
            alph(1) = alphi_values(i);
            for j=1:ngrid
                alph(2) = alphi_values(j);
                for k=1:ngrid
                    alph(3) = alphi_values(k);
                    for l=1:ngrid
                        alph(4) = alphi_values(l);
                        for m=1:ngrid
                            alph(5) = alphi_values(m);
                            for n=1:ngrid
                                alph(6) = alphi_values(n);
                                res = obj_GMM_LOMgain_univariate(alph,x,fegrid_fine,param,gx,hx,eta,e,T,ndrop,PLM,gain,p,Om,W1);
                                loss(i,j,k,l,m,n) = sum(res.^2);
                            end
                        end
                    end
                end
            end
        end
        toc
        
        filename = ['6by6loss_', todays_date];
        loss_outputs = {loss, alphi_values};
        save([filename,'.mat'],'loss_outputs')
        disp(['Saving as ' filename])
        
    end
    % Analyse the 6x6 loss
    load('6by6loss_24_Jun_2020.mat')
    loss6by6 = loss_outputs{1};
    alphi_values = loss_outputs{2};
    lossvec = loss6by6(:);
    
    load('6by6loss_pointIDs.mat') % this would take 24 min -> let the server do it too (pointIDs_server.m). It actually took
    % 0.02 sec. Ok!
    pointIDs = pointID;
    
    % 10 smallest and 10 largest losses
    [min_loss,min_idx]= mink(lossvec,10);
    [max_loss,max_idx]= maxk(lossvec,10);
    
    % Exclude the biggest 10 losses
    loss_wo_max10 = lossvec;
    loss_wo_max10(max_idx) = nan;
    
    loss_10smallest = lossvec(min_idx);
    
    % Find the alphas that cause the biggest and smallest losses
    smallestIDs = pointIDs(:,min_idx);
    largestIDs = pointIDs(:,max_idx);
    alph_vals = repmat(alphi_values,1,6,10)
    smallest_alphas = alphi_values(smallestIDs);
    largest_alphas  = alphi_values(largestIDs);
    
    % Plot full and w/o largest losses
    y = [lossvec'; loss_wo_max10'];
    figname = [this_code, '_6by6_loss_', todays_date];
    seriesnames = {'Full grid', 'W/o 10 biggest losses'};
    create_pretty_subplots(y,seriesnames,figname,print_figs)
    
    % Plot the 10 smallest
    figname = [this_code, '_6by6_loss_10smallest', todays_date];
    create_pretty_plot_x(1:length(loss_10smallest),loss_10smallest',figname,print_figs)
    
    % Plot alphas for largest and smallest losses
    figname = [this_code, '_largest_alphas', todays_date];
    legendentries = {'Largest', ''};
    create_pretty_plot_holdon(largest_alphas',legendentries,figname,print_figs)
    figname = [this_code, '_smallest_alphas', todays_date];
    legendentries = {'Smallest', ''};
    create_pretty_plot_holdon(smallest_alphas',legendentries,figname,print_figs)
    
    
end
%% 2.2) alph_true in (0,0.1) and convex
if do22==1
    command_GMM_LOMgain_univariate
end
%% 2.3) Analyze stuff from server
if do23==1
    filename = 'best_n10000_24_Jun_2020'; % simulated data, nfe = 6. full Om
    load([filename, '.mat'])
    
    true_vs_best = outputs{1};
    nsearch  = outputs{2};
    alph_opt = outputs{3};
    resnorm  = outputs{4};
    res      = outputs{5};
    Om_opt   = outputs{6};
    flag     = outputs{7}; % to check that flag seems to be always 0
    
    alph_true = true_vs_best(:,1);
    resnorm_plus = resnorm(resnorm>0);
    [min_resnorm,min_idx]= min(resnorm_plus);
    
    [alph_true, alph_opt(:,min_idx)]
    
end
%% 2.4) Add moments
if do24==1
    command_GMM_LOMgain_univariate
end

%% 2.5) Check if candidate sol on real data is robust to different starting points
if do25==1
nsearch=10;
disp(['Expected to take ', num2str(nsearch*30/60), ' minutes.'])
filename ='acf_data_11_Jun_2020'; % real data
% filename = 'acf_sim_univariate_data_21_Jun_2020'; % simulated data, nfe = 6. full Om
load([filename, '.mat'])
nfe=5;
k1min = 0;
k1max= 1;
femax = 3.5;
femin = -femax;

% Uniform random starting values
rng('default')
% b=1; a=0;
% ALPH0 = a + (b-a).*rand(nfe,nsearch);
ALPH0 = rand(nfe,nsearch);
alph_opt = ones(nfe,nsearch);
resnorm = zeros(1,nsearch);
res = zeros(45,nsearch);
Om_opt = zeros(45,nsearch);
flag = zeros(1,nsearch);

% just used to check against command_GMM_LOMgain_univariate.m
% ALPH0 =     [0.0674
%     0.0168
%     0
%     0.0168
%     0.0674]; % default*5


tic
for i=1:nsearch
    alph0 = ALPH0(:,i);
    [alph_opt(:,i), resnorm(i),res(:,i), Om_opt(:,i),flag(i)] = fun_GMM_LOMgain_univariate(acf_outputs, nfe, k1min, k1max, femin, femax, alph0);
end
toc

alph_converged = alph_opt(:,flag>0);
resnorm_converged = resnorm(flag>0);
[min_resnorm,min_idx]= min(resnorm_converged);
alph_best = alph_converged(:,min_idx);
disp('The best candidate and mean of candidates:')
[alph_best, mean(alph_converged,2)]
end
