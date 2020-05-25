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
if exist('value_outputs.mat', 'file') ~=2
    
    
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
figure
subplot(2,1,1)
plot(pgrid,pibp(1,:,1,1))
title('Pibar_t as a function of pibar_{t-1}')
subplot(2,1,2)
plot(k1grid,k1p(:,1,1,1))
title('k^{-1}_t as a function of k^{-1}_{t-1}')
sgtitle('Decision rules')

figure
subplot(2,2,1)
plot(k1grid,it(:,1,1,1))
xlabel('k^{-1}_{t-1}')
subplot(2,2,2)
plot(pgrid,it(1,:,1,1))
xlabel('pibar_{t-1}')
subplot(2,2,3)
plot(sgrid,squeeze(it(1,1,:,1)))
xlabel('r^n_t')
subplot(2,2,4)
plot(sgrid,squeeze(it(1,1,1,:)))
xlabel('u_t')
sgtitle('Policy function i_t')

close all


% suppose we have a history of states
T=100;
X = randn(4,T);
% I don't want this b/c this takes 100^4 grid
% vsim = fnval(pp,{X(1,:),X(2,:),X(3,:),X(4,:)});
% instead, I want this:
vsim = fnval(pp,X);
k1sim = X(1,:);
pibsim = X(2,:);
ssim = [X(3,:); zeros(1,T); X(4,:)];


i0 = zeros(1,T);
L = mat31_TV4(param,gx,hx,i0,pibsim,ssim,vsim);

options1.MaxFunEvals = 40000;
options1.Display = 'iter';
tv = @(i) mat31_TV4(param,gx,hx,i,pibsim,ssim,vsim);
tic
[isim, ln, flg] = fminunc(tv, i0, options1);
toc

figure
plot(isim); hold on

% another attempt at it based on Judd, 1998, p 100 (Mac).
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
% Compute this inverse in Mathematica, materials31.nb, it's called isol_{1,2}

a = [pibsim;zeros(size(pibsim));zeros(size(pibsim))];
b = gx*hx;
[fa, fb] = fafb_anal_constant_free(param,a,b,ssim,hx);
[s1, s2, s3, s4] = smat(param,hx);

s3fa = s3*fa;
s1fb = s1*fb;
s4s  = s4*ssim;
s2s  = s2*ssim;


isol1 = (kapp.^2.*sig.^2+lamx.*sig.^2).^(-1).*(kapp.^2.*s1fb.*sig+lamx.* ...
  s1fb.*sig+kapp.^2.*s2s.*sig+lamx.*s2s.*sig+kapp.*s3fa.*sig+kapp.* ...
  s4s.*sig+(-1).*((-1).*lamx.*s3fa.^2.*sig.^2+(-2).*lamx.*s3fa.* ...
  s4s.*sig.^2+(-1).*lamx.*s4s.^2.*sig.^2).^(1/2));

isol2 = (kapp.^2.*sig.^2+ ...
  lamx.*sig.^2).^(-1).*(kapp.^2.*s1fb.*sig+lamx.*s1fb.*sig+kapp.^2.* ...
  s2s.*sig+lamx.*s2s.*sig+kapp.*s3fa.*sig+kapp.*s4s.*sig+((-1).* ...
  lamx.*s3fa.^2.*sig.^2+(-2).*lamx.*s3fa.*s4s.*sig.^2+(-1).*lamx.* ...
  s4s.^2.*sig.^2).^(1/2));

plot(real(isol1))
