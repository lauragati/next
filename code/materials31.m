% materials31
% trying to implement value iteration for my model
% 22 May 2020

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

skip = 1;

%% 

[param, set, param_names, param_values_str, param_titles] = parameters_next;
sig_r = param.sig_r;
sig_u = param.sig_u;

% RE model
[fyn, fxn, fypn, fxpn] = model_NK(param);
[gx,hx]=gx_hx_alt(fyn,fxn,fypn,fxpn);

% Only begin running value iteration if optimal outputs don't exist
if exist('value_outputs.mat', 'file') ~=2
options1 = optimset('fminunc');
options1 = optimset(options1, 'TolFun', 1e-16, 'display', 'none');
% options1.UseParallel=true; % fminunc is a lot slower!
% Grids
nk = 4;
gbar = param.gbar;
k1grid = linspace(0,gbar,nk);
np = 4;
pgrid = linspace(-10,10,np);
ns = 2;
sgrid = linspace(-sig_r,sig_r,ns);
p = 0.5;
PI = [p*p, p*(1-p); (1-p)*p, (1-p)*(1-p)];

v = zeros(nk,np,ns,ns);
pp = csapi({k1grid,pgrid,sgrid,sgrid},v);
it = zeros(nk,np,ns,ns);
k1p = zeros(size(v));
pibp = zeros(size(v));
v1 = zeros(size(v));

crit=1;
iter=1;
maxiter=600;
epsi=1e-6;
datestr(now)
tic
i0 = 0; % the steady state value
while crit > epsi && iter< maxiter && crit < 1e10
    % Step 1: maximization
    for i=1:nk
        k1t_1=k1grid(i);
        for j=1:np
            pibart_1 = pgrid(j);
            for k=1:ns
                r = sgrid(k);
                for l=1:ns
                    u=sgrid(l);
                    s = [r,0,u]';
                    tv = @(i) mat31_TV3(param,gx,hx,pp,i,pibart_1,k1t_1, s,sgrid,PI);
                    it(i,j,k,l) = fminunc(tv, i0, options1);
                    % compute the value function at the maximizing i
                    [v1(i,j,k,l), pibp(i,j,k,l), k1p(i,j,k,l)] = mat31_TV3(param,gx,hx,pp,it(i,j,k,l),pibart_1,k1t_1, s,sgrid,PI);
                end
            end
        end
    end
    
    % Step 2: fitting
    pp1 = csapi({k1grid,pgrid,sgrid,sgrid},v1);
    
    % Compute stopping criterion and update
    crit = max(max(max(max(abs(v1-v)))));
    if mod(iter,10)==0 || iter==1
        disp(['Concluding iteration = ', num2str(iter)])
        disp(['Stopping criterion = ', num2str(crit)])
    end
    iter=iter+1;
    pp=pp1;
    v=v1;
end
toc
% took about 4 minutes

value_sols = {pp1,v1,it,pibp,k1p};
save('value_outputs.mat', 'value_sols')

else
    load('value_outputs.mat')
end
pp   = value_sols{1};
v    = value_sols{2};
it   = value_sols{3};
pibp = value_sols{4};
k1p  = value_sols{5};

%% 
