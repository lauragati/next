function [xn_true, xl_true, y, wl] = next_simul(param,T)
% Simulate CEMP as a function of parameters (Does the same thing as cemp_simulation, but saves output.)
pistar  = param(1);
thetbar = param(2);
gbar    = param(3);
gam     = param(4);
Gam     = param(5);
rho     = param(6);
sige    = param(7);
sigmu = param(8);
sigo1 = param(9);
sigo2 = param(10);
sigo3 = param(11);
sigo4 = param(12);
sigo5 = param(13);

k_n = 2; % number of nonlinear states
k_l = 3; % number of linear states
ne  = 2; % number of state innovations
l = 5; % number of observed variables (check that this is consistent with obs_blocksize!)

%%%------------------
%%% simulate data from CEMP's model (notes 22 June 2019)
%%%------------------

T = 500;%500;
burnin = 0;%5000;
obs_blocksize = 3; % I'm creating a dataset with no missing values

xn_true = zeros(k_n,T);
xl_true = zeros(k_l,T);
y  = zeros(l,T);
% 1.) Initialize xn0
randgbar = gbar.*rand(1,1);
k0 = randgbar.^(-1);
pibar0 = randn(1,1);
xn_t_1 = [k0; pibar0];
% evaluate initial xl0 given xn0
[~, ~, ~, fxi] = cemp_functions_matrices(param, k0, pibar0, obs_blocksize);
xi0 = fxi;
xl_t_1 = xi0;

[SIG, Sxi] = cemp_SIG_S(param); 
% generate shock sequence
wl = mvnrnd([0, 0],SIG,T+burnin)'; %  3.)  scale with chol(SIG)
e = randn(l,T+burnin); % measurement error

for t=1:T+burnin
    % 2a) Evaluate functions given xn_{t-1}
    [fk, fpibar, Apibar, fxi, Axi, Sxi, SIG, h0, hpibar, H, R] = cemp_functions_matrices(param, xn_t_1(1), xn_t_1(2), obs_blocksize);
    kt = fk;
    pibar_t = fpibar + Apibar*xl_t_1;
    xn_t = [kt,pibar_t];
    
    % Map things to SGN notation:
    [fn, An, fl, Al, Gl, h, C, Ql, R, Gn, Qn, Qln] = cemp_to_sgn_notation(fk, fpibar, Apibar, fxi, Axi, Sxi, SIG, hpibar, H, R);
    
    % 2b) Advance xl_t given xn_{t}
    xl_t = fl + Al* xl_t_1 + Gl*wl(:,t);
    
    % 2c) Evaluate yt given xl_t
    y_t = h0 + h*pibar_t + C*xl_t + chol(R)*e(:,t);
    
    if t>burnin
    tt = t-burnin;
    xn_true(:,tt) = xn_t;
    xl_true(:,tt) = xl_t;
    y(:,tt) = y_t;
    end
    % Update
    xn_t_1 = xn_t;
    xl_t_1 = xl_t;
end

