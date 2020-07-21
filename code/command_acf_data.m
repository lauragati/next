%% 1.) Get real data and filter it
    % Suppose we have real data (in logs)
    y_data = get_data;
    [ny,T] = size(y_data)
%     return
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
        lost_periods_H = h+3;

    
    % BK filter
    K=12;
    ystar_data = nan(ny,T-2*K);
    for i=1:ny
        ystar_data(i,:) = BKfilter(y_data(i,:)');
    end
        lost_periods_BK = 2*K;

    % % Plot filtered series
    create_plot_observables(c_data,{'\pi', 'x', 'i', 'E(\pi)'}, 'HP-filtered cycle', [this_code, '_HPall_', todays_date], 0)
    create_plot_observables(v_data,{'\pi', 'x', 'i', 'E(\pi)'}, 'Hamilton-filtered cycle', [this_code, '_Hamiltonall_', todays_date], 0)
    create_plot_observables(ystar_data,{'\pi', 'x', 'i', 'E(\pi)'}, 'BK-filtered cycle', [this_code, '_BKall_', todays_date], 0)


%     return
    % choose your favorite filtered data
    filt_data = ystar_data;
    lost_periods = lost_periods_BK;
    
    
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
    disp('Creating the bootstrapped autocovariances: takes about 250 sec b/c parallel')
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
    acf_outputs = {Om, Om_boot,ny,p,K, filt_data, lost_periods};
    save([filename,'.mat'],'acf_outputs')
    disp(['Saving as ' filename])