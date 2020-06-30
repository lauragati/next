% command_vfi_approx.m
% value function iteration for my model
% accelerated version, takes around 8.5 min
% now with the approximated LOM gain
% 13 June 2020

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

% Load estimation outputs
filename = 'best_n100_29_Jun_2020'; % materials35 candidate

% load the saved stuff
load([filename,'.mat'])
alph_best = output{1};
resnorm = output{2};
alph_opt = alph_best(:,1);
% grab the rest from materials35, part 2.5
nfe=5;
k1min = 0;
k1max= 1;
femax = 3.5;
femin = -femax;
% and from materials35, intro
fegrid = linspace(femin,femax,nfe);
x = cell(1,1);
x{1} = fegrid;



alph = alph_opt;



%%

[param, set, param_names, param_values_str, param_titles] = parameters_next;
sig_r = param.sig_r;
sig_u = param.sig_u;

% RE model
[fyn, fxn, fypn, fxpn] = model_NK(param);
[gx,hx]=gx_hx_alt(fyn,fxn,fypn,fxpn);

% Grids
nk = 4;
gbar = param.gbar;
% k1grid = linspace(0.0001,gbar,nk);
np = 4;
% pgrid = linspace(-0.2,0.2,np);
pgrid = linspace(-4,4,np);
ns = 2;
sgrid = linspace(-sig_r,sig_r,ns);
p = 0.5;
PI = [p*p, p*(1-p); (1-p)*p, (1-p)*(1-p)];

options1 = optimset('fminunc');
options1 = optimset(options1, 'TolFun', 1e-16, 'display', 'none');
% options1.UseParallel=true; % fminunc is a lot slower!


v = zeros(np,ns,ns,ns,ns);
pp = csapi({pgrid,sgrid,sgrid,sgrid,sgrid},v);
it = zeros(size(v));
k1p = zeros(size(v));
pibp = zeros(size(v));
v1 = zeros(size(v));

crit=1;
iter=1;
maxiter=2000; %2000
epsi=1e-6;
jj=100;
datestr(now)
tic
i0 = 0; % the steady state value
while crit > epsi && iter< maxiter && crit < 1e10
    % as VFI converges, make it start evaluating the policy function
    % more frequently
    if crit < epsi*10 && crit >= 5*epsi
        jj = 5;
    elseif crit < 2*epsi
        jj=1;
    end
    % Step 1: maximization
%     for i=1:nk
%         k1t_1=k1grid(i); % k^{-1}_{t-1} is no longer a relevant state
%         variable!
        for j=1:np
            pibart_1 = pgrid(j);
            for k=1:ns
                r = sgrid(k);
                for l=1:ns
                    u=sgrid(l);
                    s = [r,0,u]';
                    % t-1 exog states:
                    for m=1:ns
                        rt_1 = sgrid(m);
                        for n=1:ns
                            ut_1 = sgrid(n);
                            st_1 = [rt_1; 0; ut_1];
                            % only update policy every jjth iteration
                            if mod(iter,jj)==0 || iter==1
                                tv = @(i) mat31_TV3_approx(param,gx,hx,pp,i,pibart_1, s,st_1,sgrid,PI, alph,x);
                                it(j,k,l,m,n) = fminunc(tv, i0, options1);
                            end
                            % compute the value function at the maximizing i
                            [v1(j,k,l,m,n), pibp(j,k,l,m,n), k1p(j,k,l,m,n)] = mat31_TV3_approx(param,gx,hx,pp,it(j,k,l,m,n),pibart_1,s,st_1,sgrid,PI, alph,x);
                        end
                    end
                end
            end
        end
%     end
    
    % Step 2: fitting
    pp1 = csapi({pgrid,sgrid,sgrid,sgrid,sgrid},v1);
    
    % Compute stopping criterion and update
    crit = max(max(max(max(max(max(abs(v1-v)))))));
    if mod(iter,10)==0 || iter==1
        disp(['Concluding iteration = ', num2str(iter)])
        disp(['Stopping criterion = ', num2str(crit)])
    end
    iter=iter+1;
    pp=pp1;
    v=v1;
end
toc
% took about 30 minutes

value_sols = {pp1,v1,it,pibp,k1p, pgrid};
if crit < epsi
    filename = ['value_outputs_approx', nowstr, '.mat'];
    save(filename, 'value_sols')
    disp(['Saving results as ', filename])
end