% command_valfun_iter
% Attempts to find optimal i-sequence using value function iteration
% 6 April 2020
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

skip_old_stuff = 1;
%% Initialize params

[param, setp, param_names, param_values_str, param_titles] = parameters_next;
[fyn, fxn, fypn, fxpn] = model_NK(param);
[gx,hx]=gx_hx_alt(fyn,fxn,fypn,fxpn);
[ny, nx] = size(gx);
% nx is the size of s, the exogenous states. I'm going to use n_x = 4 to
% represent the number of states we're keeping track of (pibar, k^{-1}, s1(=rn), s3 (=u) )
n_x =4;

% Create grid 
N = 10; % amount of points to consider for each state variables (ngrid = N^{n_x})
[grid, ngrid] = create_grid_4states(N);

V = zeros(ngrid,1); % initial guess
tol=0.01; maxiter =300; dif = tol+1000; iter=0;
% policy function for states
policy = zeros(n_x,ngrid);
while dif > tol && iter < maxiter
    for i=1:ngrid
        % for state of the world i, the states at t take the value x
        x = grid(:,i);
        % get future states as max of current value function
%         xp = argmin stuff;
        % Evaluate tomorrow's value function as a function of the states we
        % just solved for: (get V_{t+1(x_{t+1})}, with *(-1) if we did a minimization
        Vp(i,1) = - valfun(xp);
        % Evaluate the policy function
        policy(:,i) = xp;
    end
    dif = norm(Vp-V);
    V=Vp;
    iter=iter+1;
end