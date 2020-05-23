% materials31
% trying to implement value iteration for my model

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

bet   = param.bet;
sig   = param.sig;
alph  = param.alph;
kapp  = param.kapp;
lamx  = param.lamx;
rhok  = param.rho_k;
gamk  = param.gam_k;
sig_r = param.sig_r;
sig_i = param.sig_i;
sig_u = param.sig_u;

% RE model
[fyn, fxn, fypn, fxpn] = model_NK(param);
[gx,hx]=gx_hx_alt(fyn,fxn,fypn,fxpn);
[ny, nx] = size(gx);
SIG = eye(nx).*[sig_r, sig_i, sig_u]';
eta = SIG; %just so you know

% Create convolution of parameters
b = gx*hx;
b1 = b(1,:);
b2 = b(2,:);
b3 = b(3,:);
e1 = [1,0,0]; % need to check on this!
e2 = [0,1,0];
e3 = [0,0,1];
Ibhx1 = (eye(nx)-bet*hx)^(-1);
Iabhx1 = (eye(nx)-alph*bet*hx)^(-1);
Om1 = -lamx*sig *       (sig*b1*Ibhx1 +(1-bet)*b2*Ibhx1 -sig*bet*b3*Ibhx1 +e1*sig*Ibhx1 );
Om2 = lamx*sig/(1-bet)* (sig*b1*Ibhx1 +(1-bet)*b2*Ibhx1 -sig*bet*b3*Ibhx1 +e1*sig*Ibhx1 );
Om3 = lamx *            (sig*b1*Ibhx1 +(1-bet)*b2*Ibhx1 -sig*bet*b3*Ibhx1 +e1*sig*Ibhx1 );
Om4 = kapp*sig/(1-bet) + (1-alph)*bet/(1-alph*bet);
Om5 = kapp*(sig*b1*Ibhx1 +(1-bet)*b2*Ibhx1 -sig*bet*b3*Ibhx1 +e1*sig*Ibhx1) + (1-alph)*bet*b1*Iabhx1 + kapp*alph*bet*b2*Iabhx1 + e3*Iabhx1;
Om6 = lamx*sig^2 +(kapp*sig)^2;
Om7 = lamx*sig^2/(1-bet)^2 +Om4^2;
Om8 = -lamx*sig^2/(1-bet) -kapp*sig*Om4;
Om9 = Om1 - kapp*sig*Om5;
Om10 = Om2+Om4*Om5;
Om11 = Om3 + Om5*Om5';
Om12 = Om4-1;
OM = {Om1,Om2,Om3,Om4, Om5, Om6, Om7, Om8, Om9, Om10, Om11, Om12,b1};

options1 = optimset('fminunc');
options1 = optimset(options1, 'TolFun', 1e-16, 'display', 'none');

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
maxiter=20;
epsi=1e-6;
datestr(now)
tic
i0 = 0; % the steady state value
st_1 = [0,0,0]';
while crit > epsi && iter< maxiter
    % Step 1: maximization
    for i=1:nk
        k1=k1grid(i);
        for j=1:np
            pib = pgrid(j);
            for k=1:ns
                r = sgrid(k);
                for l=1:ns
                    u=sgrid(l);
                    s = [r,0,u]';
                    tv = @(i) mat31_TV(OM,param,pp, i, pib,k1,s,st_1,sgrid,PI);
                    it(i,j,k,l) = fminunc(tv, i0, options1);
                    
                    % compute the value function at the maximizing i
                    [v1(i,j,k,l), pibp(i,j,k,l), k1p(i,j,k,l)] = mat31_TV(OM,param,pp, it(i,j,k,l), pib,k1,s,st_1,sgrid,PI);
                end
            end
        end
    end
    
    % Step 2: fitting
    pp1 = csapi({k1grid,pgrid,sgrid,sgrid},v1);
    
    % Compute stopping criterion and update
    crit = max(max(max(max(abs(v1-v)))))
    if mod(iter,50)==0
        disp(['Concluding iteration = ', num2str(iter)])
        disp(['Stopping criterion = ', num2str(crit)])
    end
    iter=iter+1
    pp=pp1;
    v=v1;
    st_1 = s;
end
toc

