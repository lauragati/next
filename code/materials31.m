% materials31
% trying to implement value iteration for my model
% 22 May 2020

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
pgrid = linspace(-10,10,np);
ns = 2;
sgrid = linspace(-sig_r,sig_r,ns);
p = 0.5;
PI = [p*p, p*(1-p); (1-p)*p, (1-p)*(1-p)];

options1 = optimset('fminunc');
options1 = optimset(options1, 'TolFun', 1e-16, 'display', 'none');
% options1.UseParallel=true; % fminunc is a lot slower!

% Only begin running value iteration if optimal outputs don't exist
if exist('value_outputs.mat', 'file') ~=2 || really_do_valiter==1
    
    v = zeros(nk,np,ns,ns,ns,ns);
    pp = csapi({k1grid,pgrid,sgrid,sgrid,sgrid,sgrid},v);
    it = zeros(size(v));
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
                        % t-1 exog states:
                        for m=1:ns
                            rt_1 = sgrid(m);
                            for n=1:ns
                                ut_1 = sgrid(n);
                                st_1 = [rt_1; 0; ut_1];
                                tv = @(i) mat31_TV3(param,gx,hx,pp,i,pibart_1,k1t_1, s,st_1,sgrid,PI);
                                it(i,j,k,l,m,n) = fminunc(tv, i0, options1);
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
    % took about 4 minutes
    
    value_sols = {pp1,v1,it,pibp,k1p};
    if crit < epsi
        save('value_outputs.mat', 'value_sols')
    end
else
%     load('value_outputs.mat')
load('value_outputs_server.mat')
end
pp   = value_sols{1};
v    = value_sols{2};
it   = value_sols{3};
pibp = value_sols{4};
k1p  = value_sols{5};

%%
% suppose we have a history of states from parametric expectations
load('inputs.mat')
seriesnames = {'\pi', 'x','i'};
invgain = {'Inverse gain'};

e = output{1};
ysim7 = output{2};
k7    = output{3};
phi7  = output{4};
seq_opt = output{5};
i_opt = seq_opt(3,:);

% take the states from the results of parametric E
T=length(e)-2; % drop the first and last obs
k1sim = 1./k7(2:end-1);
pibsim = squeeze(phi7(1,1,2:end-1))';
ssim = e(:,2:end-1);
ssimt_1 = e(:,1:end-2);
X = [k1sim; pibsim; ssim(1,:); ssim(3,:); ssimt_1(1,:); ssimt_1(3,:)];

% I don't want this b/c this takes 100^4 grid
% vsim = fnval(pp,{X(1,:),X(2,:),X(3,:),X(4,:)});
% instead, I want this:
vsim = fnval(pp,X);
vsim = 1000*ones(size(vsim)); % replace as long as vsim isn't accurate

% attempt to back out policy from the converged V based on Judd, 1998, p 100 (Mac).
bet = param.bet;
lamx = param.lamx;
sig  = param.sig;
kapp = param.kapp;
% compute the period loss from time-invariant policy and value function
Pu = (1-bet)*vsim; % Judd (12.3.5.)
% now I'd need to invert the period loss function for pi, x, i
% Note that Pu = pi^2 + lamx*x^2, and
% x = -sig.*i + s1*fb + s2*s;
% pi = kapp*x + s3*fa + s4*s;, so

a = [pibsim;zeros(size(pibsim));zeros(size(pibsim))];
b = gx*hx;
[fa, fb] = fafb_anal_constant_free(param,a,b,ssim,hx);
[s1, s2, s3, s4] = smat(param,hx);

s3fa = s3*fa;
s1fb = s1*fb;
s4s  = s4*ssim;
s2s  = s2*ssim;

L = Pu;
% Compute this inverse in Mathematica, materials31.nb, it's called isol_{1,2}
isol_mma = (kapp.^2.*sig.^2+lamx.*sig.^2).^(-1).*(kapp.^2.*s1fb.*sig+lamx.* ...
    s1fb.*sig+kapp.^2.*s2s.*sig+lamx.*s2s.*sig+kapp.*s3fa.*sig+kapp.* ...
    s4s.*sig+(-1).*(kapp.^2.*L.*sig.^2+L.*lamx.*sig.^2+(-1).*lamx.* ...
    s3fa.^2.*sig.^2+(-2).*lamx.*s3fa.*s4s.*sig.^2+(-1).*lamx.*s4s.^2.* ...
    sig.^2).^(1/2));

% or compute it numerically here
options1.MaxFunEvals = 40000;
options1.Display = 'iter';
options1.TolFun = 1e-20;
options1.TolX = 1e-20;


i0 = i_opt;
l = @(i) period_loss_i(param,hx,i,fa,fb,ssim,L);
tic
[isol_matlab, ln, flg] = fminunc(l, i0, options1);
toc

% solve it initializing away from PE sol
l = @(i) period_loss_i(param,hx,i,fa,fb,ssim,L);
tic
[isol_matlab_initrand, lnr, flgr] = fminunc(l, i0+100*rand(1,T), options1);
toc

figure
plot(real(isol_mma));  hold on
plot(isol_matlab)
plot(isol_matlab_initrand)
plot(i_opt, 'k--')
legend('mma', 'matlab, init at PE','matlab, init rand', 'sol PE')
title('')

% figure
% plot(i_opt, 'k'); hold on
% plot(isol_matlab, '--')
% legend('sol from parametric E', 'sol from VFI')
% title('bingo')

%create_plot_observables([i_opt; isol_matlab],{'PE','VFI',}, 'Optimal i-sequence given a simulation, initalized at i^{PE} ', 'pe_vfi_initPE', print_figs)
%create_plot_observables([i_opt; isol_matlab_initrand],{'PE','VFI',}, 'Optimal i-sequence given a simulation, initalized at i^{PE}+10*rand ', 'pe_vfi_initrand', print_figs)

