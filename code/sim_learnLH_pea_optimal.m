function [xsim, ysim] = sim_learnLH_pea_optimal(alph,x,param,gx,hx,eta, PLM, gain,T,ndrop,e,g_fe,knowTR,mpshock)
% simulate model given shocks and solving for the PEA-optimal interest rate
% sequence
% With outputs of estimation of univariate approximated LOM gain / anchoring function
% takes around 5*N min
% 30 July 2020



%% Initial sequence

%Select exogenous inputs
s_inputs = [1;1;1]; % pi, x, i

% find indeces and number of input sequences
i_inputs = find(s_inputs); % index of inputs series in y
n_inputs = sum(s_inputs); % the number of input series

v=zeros(4,T); % need to deal with the darn measurement error later

% an initial simulation using the Taylor rule
[~, y0, ~, ~, ~, ~, FEt_10, ~] = sim_learnLH_clean_approx_univariate(alph,x,param,gx,hx,eta, PLM, gain, T+ndrop,ndrop,e,v, 1,mpshock);


% Note: I'm not inputting anything exogenous for period t=1 b/c that
% just causes errors that by construction fsolve can't close
seq0 = zeros(n_inputs,T);
seq0(:,2:end) = y0(i_inputs,2:end);


%% Parameterized expectations

% %Optimization Parameters
options = optimoptions('fsolve', 'TolFun', 1e-9, 'display', 'iter');%, 'MaxFunEvals', 4000);
options.UseParallel=true;
% initDamping = initial value of Levenberg-Marquardt lambda.

% 7.) Parameterized expectations approach
% initialize at Taylor rule
seq0crop = [seq0(:,2:end-1);FEt_10(1,2:end-1)]; % input jumps and Fe(pi)
% Or detach from Taylor rule
seq0crop = rand(size(seq0crop));
% initialize beta-coefficients
bet0 = zeros(4,3); % 4 states x 3 powers % INITIALIZE AT ZEROS INSTEAD OF ONES AND IT FINDS SOL
bet = bet0(:);

% return
fegrid = x{1};

% Evaluate residuals once
resids = objective_seq_clean_parametricE_approx(seq0crop,bet,n_inputs,param,gx,hx,eta,PLM,gain,T,ndrop,e, alph, x,fegrid, g_fe, knowTR);
% disp('Initial residuals IS, PC, TC and A7')
% disp(num2str(resids))

% return

dbstop if warning

maxiter=12;
BET = zeros(length(bet),maxiter);
iter=0;
crit=1;
start = now;
datestr(now)
while crit > 1e-6 && iter < maxiter
    iter=iter+1
    BET(:,iter) = bet; % storing betas
    % Now solve model equations given conjectured E 
    objh = @(seq) objective_seq_clean_parametricE_approx(seq,bet,n_inputs,param,gx,hx,eta,PLM,gain,T,ndrop,e, alph,x,fegrid, g_fe, knowTR);
    
    tic
    [seq_opt, ~, flag] = fsolve(objh,seq0crop, options);
    toc
    % seq_opt-[seq0(:,2:end-1);FEt_10(1,2:end-1)]
    if flag==1 % If fsolve converged to a root
    % Projection step: Recover v, compute analogues E, update beta
    bet1 = projection_approx(seq_opt,n_inputs,param,gx,hx,eta,PLM,gain,T,ndrop,e,alph,x,fegrid, g_fe, knowTR);
    crit=max(abs(bet-bet1));
    bet=bet1;
    seq0crop = seq_opt; % try to accelerate 
    else
        disp('No sol, stopping PEA.')
        return
    end
end
endt = now;
elapsed_seconds = etime(datevec(endt), datevec(start));
disp(['Elapsed: ' num2str(elapsed_seconds), ' sec.'])

[xsim, ysim] = sim_learnLH_clean_given_seq3_approx(param,gx,hx,eta,PLM, gain, T,ndrop,e,seq_opt,n_inputs, alph,x,fegrid, g_fe, knowTR);

disp('Done.')