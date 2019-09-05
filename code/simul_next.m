function [xn_true, xl_true, y, e] = simul_next(param,set,T, burnin)
% Simulate the CEMP-Preston mixed model. Builds on cemp_simul_fun.m.
% 4 Sep 2019

bet  = param(1);
sig  = param(2);
alph = param(3);
kapp = param(4);
psi_x  = param(5);
psi_pi = param(6);
w    = param(7);
gbar = param(8);
thetbar = param(9);
rho_r = param(10);
rho_i = param(11);
rho_u = param(12);
sig_r = param(13);
sig_i = param(14); 
sig_u = param(15); 
ne = set(16); % number of state innovations
nxnl = set(17);
nxl = set(18);
ny = set(19);
nxnl2 = set(20);
P = eye(ne).*[rho_r, rho_i, rho_u]';
SIG = eye(ne).*[sig_r, sig_i, sig_u]';

k_n = nxnl; % number of nonlinear states
k_l = nxl; % number of linear states

[A1, A2, A3] = matrices_next(param, set);

xn_true = zeros(k_n,nxnl2,T);
xl_true = zeros(k_l,T);
y  = zeros(ny,T);
% 1.) Initialize xn0
randgbar = gbar.*rand(nxnl,1);
k0 = randgbar.^(-1);
zbar0 = randn(nxl,1);
s0 = zeros(ne,1);
kt_1 = k0;
zbart_1 = zbar0;
st_1 = s0;

% generate shock sequence
e = mvnrnd([0, 0, 0],SIG,T+burnin)'; 

for t=1:T+burnin
    % Exog shock process
    if t==1
        st = e(:,t);
    else
        st = P*st_1 + e(:,t);
    end
    % 2a) Evaluate functions given xn_{t-1}
    fk = functions_next(param,set,zbart_1,st_1, kt_1);
    kt = fk;
    zbart = (eye(ne) + kt.^(-1).*(A1/(1-alph*bet) +A2/(1-bet)-eye(ne)) ) * zbart_1 ...
        + kt.^(-1).*(A1*(eye(ne)-alph*bet*P)^(-1) +A2*(eye(ne)-bet*P)^(-1) +A3 - eye(ne)) * st;
    xn_t = [1./kt,zbart]; % save inverse gains
    
    % 2b) Advance xl_t given xn_{t}
    zt = (A1/(1-alph*bet) +A2/(1-bet)-eye(ne))*zbart + (A1*(eye(ne)-alph*bet*P)^(-1) +A2*(eye(ne)-bet*P)^(-1) +A3)*st;
    xl_t = zt;
    
    % 2c) Evaluate yt given xl_t - right now the observables are just zt
    y_t = zt;
    
    if t>burnin
    tt = t-burnin;
    xn_true(:,:,tt) = xn_t;
    xl_true(:,tt) = xl_t;
    y(:,tt) = y_t;
    end
    % Update
    kt_1 = kt;
    zbart_1 = zbart;
    st_1 = st;
end

