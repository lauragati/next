% materials33
% begin to estimate anchoring function
% return to the old GMM code, but write a new objective function
% 5 June 2020

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
really_do_valiter = 0;

%% Takes it over from command_GMM_anchoring_function

% command_GMM_anchoring_function.m
% a GMM implementation of estimating the anchoring function g
% created in materials22.m and replaces that file.
% 26 March 2020

clearvars
close all
clc

% Add all the relevant paths and grab the codename
this_code = mfilename;
[current_dir, basepath, BC_researchpath,toolpath,export_figpath,figpath,tablepath,datapath,tryouts_path,inputsRyan_path] = add_paths;
todays_date = strrep(datestr(today), '-','_');

% Variable stuff ---
print_figs        = 0;
stop_before_plots = 0;
skip_old_plots    = 0;
output_table = print_figs;

redo_data_load_and_bootstrap = 0;

if redo_data_load_and_bootstrap==1 % takes about 140 sec
    %% 1.) Get real data and filter it
    % Suppose we have real data (in logs)
    y_data = get_data;
    [ny,T] = size(y_data)
    
    % Filter the data
    % HP filter
    g_data = nan(size(y_data));
    c_data = nan(size(y_data));
    for i=1:ny
        [g_data(i,:),c_data(i,:)] = HPfilter(y_data(i,:)');
    end
    
    % Hamilton filter
    h=8;
    v_data = nan(ny,T-4-h+1);
    for i=1:ny
        [v_data(i,:)] = Hamiltonfilter(y_data(i,:)');
    end
    
    
    % BK filter
    K=12;
    ystar_data = nan(ny,T-2*K);
    for i=1:ny
        ystar_data(i,:) = BKfilter(y_data(i,:)');
    end
    
    % % Plot filtered inflation
    create_plot(1:numel(c_data(1,:)),c_data(1,:),'HP-filtered cycle',[this_code, '_HP'],print_figs,'PCE Inflation')
    create_plot(1:numel(v_data(1,:)),v_data(1,:),'Hamilton-filtered cycle',[this_code, '_Hamilton'],print_figs,'PCE Inflation')
    create_plot(1:numel(ystar_data(1,:)),ystar_data(1,:),'BK-filtered cycle',[this_code, '_BK'],print_figs,'PCE Inflation')
    
    
    % choose your favorite filtered data
    filt_data = ystar_data;
    
    
    
    % compute moments, Om, as the autocovariances of the data for lags
    % 0,1,...,K
    K=4;
    % Take the initial data, estimate a VAR
    max_lags = 16;
    [AIC,BIC,HQ] = aic_bic_hq(filt_data',max_lags);
    p =min([AIC,BIC,HQ]);
    % A is the impact matrix, identified via Cholesky, B is the beta_ols, res are
    % the residuals, sigma is the estimated VC matrix.
    [A,B,res,sigma] = sr_var(filt_data', p);
    
    % Rewrite the VAR(p) as VAR(1) (see Hamilton, p. 273, Mac)
    c = B(1,:); % coefficients of the constant
    PHI = B(2:end,:)';
    F = [PHI; [eye(ny*(p-1)), zeros(ny*(p-1),ny)]];
    Q = [[sigma, zeros(ny,ny*(p-1))]; zeros(ny*(p-1),ny*p)];
    % check sizes
    np = ny*p;
    sizeF = size(F) == [np,np];
    sizeQ = size(Q) == [np,np];
    
    vecSig = (eye(np^2)-kron(F,F))\vec(Q);
    % VC matrix of data y
    Gamj = zeros(ny,ny,K+1);
    Sig = reshape(vecSig,np,np);
    for j=0:K
        % jth Autocov of data y, still as a VAR(1)
        Sigj = F^j * Sig;
        % jth Autocov of data y, finally as a VAR(p)
        Gamj(:,:,j+1) = Sigj(1:ny,1:ny);
        % This only takes the covariances wrt the same var: (e.g Cov(x_t,x_{t-1}))
    end
    % moments vector
    Om = vec(Gamj);
    
    
    
    %% 2.) Bootstrap the data, and get variance of moments (autocovariances from 0 to lag K) from the bootstrapped sample, takes more than 10 min!
    % This section is inspired by main_file_SVAR_just_IT_controllingNEWS.m in my work with Marco
    
    % Resample the residuals and use beta_ols from VAR to create nboot new samples.
    nboot =10000;
    q=16; % blocksize
    nburn = 1000;
    which_correction ='blocks';
    disp(datestr(now))
    disp('Creating the bootstrapped sample: takes about 50 sec')
    tic
    [dataset_boot] = data_boot(B,nburn,res,nboot,which_correction,q);
    toc
    % Autocov matrix from bootstrapped sample for lags 0,...,K
    K = 4;
    Gamj = zeros(ny,ny,K+1,nboot);
    Om_boot = zeros(length(Om),nboot);
    tic
    disp(datestr(now))
    disp('Creating the bootstrapped autocovariances: takes about 30 sec b/c parallel')
    parfor i=1:nboot
        Gamj_booti = zeros(ny,ny,K+1);
        %         max_lags = 16;
        %         [AIC,BIC,HQ] = aic_bic_hq(squeeze(dataset_boot(:,:,i)),max_lags);
        %         p = min([AIC,BIC,HQ]); % lag selection (p) is the lag
        [A,B,res,sigma] = sr_var(squeeze(dataset_boot(:,:,i)), p);
        
        % Rewrite the VAR(p) as VAR(1) (see Hamilton, p. 273, Mac)
        c = B(1,:); % coefficients of the constant
        PHI = B(2:end,:)';
        F = [PHI; [eye(ny*(p-1)), zeros(ny*(p-1),ny)]];
        v = [res'; zeros(ny*(p-1),size(res',2))];
        Q = [[sigma, zeros(ny,ny*(p-1))]; zeros(ny*(p-1),ny*p)];
        % check sizes
        np = ny*p;
        sizeF = size(F) == [np,np];
        sizeQ = size(Q) == [np,np];
        
        vecSig = (eye(np^2)-kron(F,F))\vec(Q);
        % VC matrix of data y
        Sig = reshape(vecSig,np,np);
        for j=0:K
            % jth Autocov of data y, still as a VAR(1)
            Sigj = F^j * Sig;
            % jth Autocov of data y, finally as a VAR(p)
            Gamj_booti(:,:,j+1) = Sigj(1:ny,1:ny);
        end
        % gather the ACF of each bootstrapped sample
        Gamj(:,:,:,i) = Gamj_booti;
        Om_boot(:,i) = vec(Gamj_booti);
    end
    toc
    filename = ['acf_data_', todays_date];
    acf_outputs = {Om, Om_boot,ny,p,K, filt_data};
    save([filename,'.mat'],'acf_outputs')
    
end % end of doing the whole data loading and bootstrapping again

%% Compute weighting matrix and do GMM
filename ='acf_data_09_Jun_2020';
load([filename, '.mat'])
Om = acf_outputs{1}; % this is the moments vector
Om_boot = acf_outputs{2}; % moments vectors in bootstrapped samples
ny = acf_outputs{3};
p = acf_outputs{4};
K = acf_outputs{5};
filt_data = acf_outputs{6};
T = length(filt_data);

% weighting matrix for GMM
% note: 2nd argument of var(X,W,DIM) signifies whether to normalize by N-1
% or N. 0 -> N-1, 1 -> N. The third argument, DIM, says along which
% dimension to take the variance.
W = diag(var(Om_boot,0,2));
W1 = W^(-1);

[param, set, param_names, param_values_str, param_titles] = parameters_next;
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

% Model selection
%%%%%%%%%%%%%%%%%%%
PLM = constant_only;
gain = again_critsmooth;
%%%%%%%%%%%%%%%%%%%

% Size of cross-section
% we're not doing a whole cross-section here
ndrop = 100; % 100

% gen all the N sequences of shocks at once.
rng(0)
e = randn(ne,T+ndrop);

% Fmincon
%Optimization Parameters
options = optimset('fmincon');
options = optimset(options, 'TolFun', 1e-9, 'display', 'iter');

% Do an initial approx of the anchoring function to initialize the coeffs
ng = 10;
% grids for k^(-1)_{t-1} and f_{t|t-1}
k1grid = linspace(0.001,param.gbar,ng);
fegrid = linspace(-5,5,ng);
% values for k^{-1}_t for the grid
k = zeros(ng,ng);
for i=1:ng
    for j=1:ng
        k(i,j) = fk_smooth_pi_only(param,fegrid(j), 1./k1grid(i));
    end
end
k1 = 1./k;
% map to ndim_simplex
x = cell(2,1);
x{1} = k1grid;
x{2} = fegrid;
[xxgrid, yygrid] = meshgrid(k1grid,fegrid);
alph = ndim_simplex(x,[xxgrid(:)';yygrid(:)'],k1);
% just a check for use
xx = [k1grid(1); fegrid(1)];
out = ndim_simplex_eval(x,xx,alph);
if out==k1(1)
    disp('bingo')
end

alph0 = alph;
ub = alph0+0.1;
lb = alph0-0.1;
% %Compute the objective function one time with some values
loss = obj_GMM_LOMgain(alph0,x,param,e,T,ndrop,PLM,gain,Om,W1);


tic
%Declare a function handle for optimization problem
objh = @(alph) obj_GMM_LOMgain(alph,x,param,e,T,ndrop,PLM,gain,Om,W1);
[alph_opt, loss_opt] = fmincon(objh, alph0, [],[],[],[],lb,ub,[],options);
% fmincon(FUN,X0,A,B,Aeq,Beq,LB,UB,NONLCON,OPTIONS)
toc

alph0-alph_opt

