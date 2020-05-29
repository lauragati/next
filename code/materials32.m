% materials32
% trying to map value iteration results to parametric expectations results
% trying to accelerate vfi (only eval every jjth policy function) --> 50
% min!
% 28 May 2020

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
k1grid = linspace(0,gbar,nk);
np = 4;
pgrid = linspace(-0.2,0.2,np);
ns = 2;
sgrid = linspace(-sig_r,sig_r,ns);
p = 0.5;
PI = [p*p, p*(1-p); (1-p)*p, (1-p)*(1-p)];

options1 = optimset('fminunc');
options1 = optimset(options1, 'TolFun', 1e-16, 'display', 'none');
% options1.UseParallel=true; % fminunc is a lot slower!

% Only begin running value iteration if optimal outputs don't exist
if exist('value_outputs.mat', 'file') ~=2 && exist('value_outputs_server.mat', 'file') ~=2 || really_do_valiter==1
    
    v = zeros(nk,np,ns,ns,ns,ns);
    pp = csapi({k1grid,pgrid,sgrid,sgrid,sgrid,sgrid},v);
    it = zeros(size(v));
    k1p = zeros(size(v));
    pibp = zeros(size(v));
    v1 = zeros(size(v));
    
    crit=1;
    iter=1;
    maxiter=2000;
    epsi=1e-3;
    jj=100;
    datestr(now)
    tic
    i0 = 0; % the steady state value
    while crit > epsi && iter< maxiter && crit < 1e10
        % as VFI converges, make it start evaluating the policy function
        % more frequently
        if crit < epsi*10 && crit >= 5*epsi
            jj = 5;
        elseif crit < 5*epsi
            jj=1;
        end
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
                        % t-1 exog states:
                        for m=1:ns
                            rt_1 = sgrid(m);
                            for n=1:ns
                                ut_1 = sgrid(n);
                                st_1 = [rt_1; 0; ut_1];
                                % only update policy every jjth iteration
                                if mod(iter,jj)==0 || iter==1
                                    tv = @(i) mat31_TV3(param,gx,hx,pp,i,pibart_1,k1t_1, s,st_1,sgrid,PI);
                                    it(i,j,k,l,m,n) = fminunc(tv, i0, options1);
                                end
                                % compute the value function at the maximizing i
                                [v1(i,j,k,l,m,n), pibp(i,j,k,l,m,n), k1p(i,j,k,l,m,n)] = mat31_TV3(param,gx,hx,pp,it(i,j,k,l),pibart_1,k1t_1,s,st_1,sgrid,PI);
                            end
                        end
                    end
                end
            end
        end
        
        % Step 2: fitting
        pp1 = csapi({k1grid,pgrid,sgrid,sgrid,sgrid,sgrid},v1);
        
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
    % took about 50 minutes
    
    value_sols = {pp1,v1,it,pibp,k1p};
    if crit < epsi
        save('value_outputs.mat', 'value_sols')
        disp('saving results...')
    end
else
%     load('value_outputs.mat')
    load('value_outputs_server.mat')
%     load('value_outputs_server32_accelerated.mat')
end


pp   = value_sols{1};
v    = value_sols{2};
it   = value_sols{3};
pibp = value_sols{4};
k1p  = value_sols{5};

%%
% take the history of states from parametric expectations
load('inputs.mat')
e = output{1};
ysim7 = output{2};
k7    = output{3};
phi7  = output{4};
seq_opt = output{5};
i_pe = seq_opt(3,:);

% compile state vector
T=length(e)-2; % drop the first and last obs
k1sim = 1./k7(2:end-1);
pibsim = squeeze(phi7(1,1,2:end-1))';
ssim = e(:,2:end-1);
ssimt_1 = e(:,1:end-2);
X = [k1sim; pibsim; ssim(1,:); ssim(3,:); ssimt_1(1,:); ssimt_1(3,:)];

% approximate policy function
ppi = csapi({k1grid,pgrid,sgrid,sgrid,sgrid,sgrid},it);
i_vi = fnval(ppi,X);

policies = [i_pe; i_vi];
create_simple_plot(policies,{'PE', 'VFI'},'Policy: i(X)',[this_code, 'accelerated_policy'],print_figs)

working_figure_name_was = [this_code, '_policy']; % DONT USE THIS!!!