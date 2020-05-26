% materials31_server.m
% value iteration for my model on the server
% do it and save results
% 26 May 2020

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
    save('value_outputs_server.mat', 'value_sols')
end