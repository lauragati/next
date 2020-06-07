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
[current_dir, basepath, BC_researchpath,toolpath,export_figpath,figpath,tablepath,datapath] = add_paths;
todays_date = strrep(datestr(today), '-','_');

% Variable stuff ---
print_figs        = 0;
stop_before_plots = 0;
skip_old_plots    = 0;
output_table = print_figs;

redo_data_load_and_bootstrap = 0;

if redo_data_load_and_bootstrap==1
    %% 1.) Get real data and filter it
    % Suppose we have real data (in logs)
    y_data = get_data;
    [ny,T] = size(y_data)
    
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
    PLM = constant_only_pi_only;
    gain = again_critsmooth;
    %%%%%%%%%%%%%%%%%%%
    
    % Size of cross-section
    N = 1; % we're not doing a whole cross-section here
    burnin = 100; % 100
    
    % gen all the N sequences of shocks at once.
    rng(0)
    eN = randn(ne,T+burnin,N);
    
    
    % Filter both using the same filter, at this point only the data
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
    
    %% 2.) Bootstrap the data, and get autocovariances from 0 to lag K from the bootstrapped sample, takes more than 10 min!
    %%%%%%%%%%%%%%% 5 JUNE ONWARDS: ESTIMATE A VAR ON THE FILTERED DATA AND %%%%%%%%%%%%%%%
    %%%%%%%%%%%%%%% BOOTSTRAP THE VAR RESIDUALS                             %%%%%%%%%%%%%%%
    % This section is inspired by main_file_SVAR_just_IT_controllingNEWS.m in
    % my work with Marco
    
    % Take the initial data, estimate a VAR and bootstrap it
    max_lags = 16;
    [AIC,BIC,HQ] = aic_bic_hq(filt_data',max_lags);
    nlags =min([AIC,BIC,HQ]);
    % A is the impact matrix, identified via Cholesky, B is the beta_ols, res are
    % the residuals, sigma is the estimated VC matrix.
    [A,B,res,sigma] = sr_var(filt_data', nlags);
    
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
    tic
    disp(datestr(now))
    disp('Creating the bootstrapped autocovariances: takes about 10 min')
    for i=1:nboot
        max_lags = 16;
        [AIC,BIC,HQ] = aic_bic_hq(squeeze(dataset_boot(:,:,i)),max_lags);
        p = min([AIC,BIC,HQ]); % lag selection (p) is the lag
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
            Gamj(:,:,j+1,i) = Sigj(1:ny,1:ny);
            % This only takes the covariances wrt the same var: (e.g Cov(x_t,x_{t-1}))
            acf_proper(:, j+1,i) = diag(Gamj(:,:,j+1,i));
            % Note: using the full autocov matrix gives ny*ny*(K+1) moments (45),
            % while using only the "idiosyncratic" or proper autocovs gives you
            % ny*(K+1) (here 45). So you're more overidentified if you use Gamj,
            % while less computationally heavy if you use acf_data_alt.
        end
    end
    toc
    meanGamj = mean(Gamj,4);
    mean_acf_proper = mean(acf_proper,3);
    filename = ['acf_data_', todays_date];
    acf_outputs = {meanGamj, mean_acf_proper,ny,p,K, filt_data};
    save([filename,'.mat'],'acf_outputs')
    
end % end of doing the whole data loading and bootstrapping again

filename ='acf_data_07_Jun_2020';
load([filename, '.mat'])
Gamj = acf_outputs{1}; % this is the full ACF matrix
acf_proper = acf_outputs{2}; % this is just the ACF of the individual variables on themselves only (it's a perfect subset of Gamj - I've checked, they agree!)
ny = acf_outputs{3};
p = acf_outputs{4};
K = acf_outputs{5};
filt_data = acf_outputs{6};

% should prefer the full ACF matrix
W = diag(vec(Gamj));
W1 = W^(-1);

% over the "proper" subset
W_proper = diag(vec(acf_proper));
W1_proper = W_proper^(-1);


%  CONT HERE WITH NEW FMINCON OBJECTIVE FUNCTION

% % % 6.) Do fmincon -- OLD STUFF, % try estimating d and c
% % %Optimization Parameters
% % options = optimset('fmincon');
% % options = optimset(options, 'TolFun', 1e-9, 'display', 'iter');
% %
% % varp0 = [10,10];
% % ub = [100,100];
% % lb = [0.1,0.1];
% % % %Compute the objective function one time with some values
% % % loss = objective_ACstructure(varp0,set,eN,T,N,burnin,PLM,gain, vec(acf_data),W1);
% %
% % tic
% % %Declare a function handle for optimization problem
% % objh = @(varp) objective_ACstructure(varp,set,eN,T,N,burnin,PLM,gain, vec(acf_data));
% % [par_opt, loss_opt] = fmincon(objh, varp0, [],[],[],[],lb,ub,[],options);
% % % fmincon(FUN,X0,A,B,Aeq,Beq,LB,UB,NONLCON,OPTIONS)
% % toc
% %
% % par_opt
% %
% % loss_opt
% %
% % % 7). Simulate the model using the optimal d:
% % param.d = par_opt(1);
% % param.c = par_opt(2);
% % [y,k] = fun_sim_anchoring(param,T,N, burnin,eN,PLM,gain);
% %
% % if size(filt_data,2) == size(c_data,2)
% %     create_plot((1:numel(k)),1./k','k^{-1}',[this_code, '_gain_par_opt_HP', '_', todays_date],print_figs,'Inverse gain given estimated params')
% % elseif size(filt_data,2) == size(v_data,2)
% %     create_plot((1:numel(k)),1./k','k^{-1}',[this_code, '_gain_par_opt_Hamilton', '_', todays_date],print_figs,'Inverse gain given estimated params')
% % elseif size(filt_data,2) == size(ystar_data,2)
% %     create_plot((1:numel(k)),1./k','k^{-1}',[this_code, '_gain_par_opt_BK', '_', todays_date],print_figs,'Inverse gain given estimated params')
% % end