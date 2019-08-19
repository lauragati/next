% command_preston3
% Compare dynamics of Preston's system (eq 18-19) to Bullard & Mitra's (eq 13-14)
% 8 August 2019
clearvars

% Add all the relevant paths
current_dir = pwd;
cd ../ % go up 1 levels
basepath = pwd;
cd .. % go up another level to BC_Research
BC_researchpath = pwd;
toolpath = [BC_researchpath '/matlab_toolbox'];
export_figpath = [toolpath '/Export_Fig'];
figpath = [basepath '/figures'];
tablepath = [basepath '/tables'];
datapath = [basepath '/data'];
lorenzonipath = [current_dir '/lorenzoni_baseline'];
ryan_learnpath = [current_dir, '/Learning_models_from_Ryan/Learning_Example_Money'];

cd(current_dir)

addpath(basepath)
addpath(toolpath)
addpath(export_figpath)
addpath(figpath)
addpath(datapath)
addpath(lorenzonipath)
addpath(ryan_learnpath)
% Variable stuff ---
print_figs    = 0;
do_old_plots  = 0;
if print_figs ==1
    output_table  = 1;
else
    output_table =0;
end
fs=20; % fontsize

clc

%% Starting values for all four models
T = 2000;
param = param_preston;
n = param.n;
f       = param.f;
C       = param.C; % IB is on top, RN is below
alph    = param.alph;
bet     = param.bet;
sig     = param.sig;
kapp    = param.kapp;
psi_x   = param.psi_x;
psi_pi  = param.psi_pi;
w_small = param.w_small;

% Generate shock sequence
rng(0)
eta = eye(n);
eps_s = randn(n,T);
s = zeros(n,T);
s(:,1) = eps_s(:,1);
for t=2:T
    s(:,t) = C*s(:,t-1) + eta*eps_s(:,t);
end
% s = randn(n,T);

% solve RE version of model to initialize
[fyn, fxn, fypn, fxpn] = model_preston(param);
[gx,hx]=gx_hx_alt(fyn,fxn,fypn,fxpn);
[ny, nx] = size(gx);
l = ny; % l is the number of regressands
k = n+1; % k is the number of regressors
gl = [zeros(ny,1) gx*hx];
hl = [zeros(nx,1), hx];
[~,sigx] = mom(gx,hx,eta*eta');
R = eye(nx+1); R(2:end,2:end) = sigx;
phi0 = gl';
R0 = R;

% learning convergence tolerance
tol= 1e-4;

%% Try Ryan's thing
[x2,y2] = sim_learn(gx,hx,fxpn,fxn,fypn,fyn,eta,T,0,zeros(nx,1),0.1);

%% Bullard & Mitra economy
[BM1, BM2] = bullard_mitra_matrices(param);

tic
[z_BM, learnin_BM, diff_BM, a_BM, b_BM] = sim_euler_eq_learning(param,gx,R0,phi0,BM1,BM2,s,tol);
toc

%%  Evolution of economy, Preston
[A1_Preston, A2_Preston, A3_Preston] = preston_matrices(param);

tic
[z_Preston, learnin_Preston, diff_Preston, a_Preston, b_Preston] = sim_long_horizon_learning(param,R0,phi0,A1_Preston, A2_Preston, A3_Preston,s,tol);
toc

%% Preston, Peter's ALM
[A1_Peter, A2_Peter, A3_Peter] = preston_peters_matrices(param);
tic
[z_Peter, learnin_Peter, diff_Peter, a_Peter, b_Peter, FA, FB, FY] = sim_long_horizon_learning(param,R0,phi0,A1_Peter, A2_Peter, A3_Peter,s,tol);
toc

% Peter's system, equation by equation
tic
[zz] = sim_long_horizon_learning2(param,R0,phi0,s,tol); % is the same as Peter's to 10e-11
toc

%% Check if simulations are correct:

Ez = zeros(l,T);
for t=1:T
    if t==1
        at_1 = phi0(1,:)';
        bt_1 = phi0(2:end,:)';
        Ez(:,t) = at_1 + bt_1*s(:,t);
    else
        Ez(:,t) = a_BM(:,t-1) + squeeze(b_BM(:,:,t-1))*s(:,t);
    end
end
Ex = Ez(2,:); Epi = Ez(1,:);
% Bullard Mitra
eq13 = -z_BM(2,:) + Ex - sig*(z_BM(3,:) - Epi) + s(2,:);
eq14 = -z_BM(1,:) + kapp*z_BM(2,:) + bet*Epi;
eq27 = -z_BM(3,:) + psi_pi*z_BM(1,:) + psi_x*z_BM(2,:) + s(1,:);

% Preston
X1_Preston = zeros(1,T);
X2_Preston = zeros(1,T);
for t=1:T
    if t==1
        at_1 = phi0(1,:)';
        bt_1 = phi0(2:end,:)';
        fb = at_1/(1-bet) + bt_1*(eye(n)-bet*C)^(-1)*s(:,t);
        fa = at_1/(1-alph*bet) + bt_1*(eye(n)-alph*bet*C)^(-1)*s(:,t);
    else
        fb = a_Preston(:,t-1)/(1-bet) + squeeze(b_Preston(:,:,t-1))*(eye(n)-bet*C)^(-1)*s(:,t);
        fa = a_Preston(:,t-1)/(1-alph*bet) + squeeze(b_Preston(:,:,t-1))*(eye(n)-alph*bet*C)^(-1)*s(:,t);
    end
    fy = f'*(eye(n)-bet*C)^(-1)*s(:,t);
    
    % these are the LR expectations on RHS of 18 and 19:
    X1_Preston(t) = (1-bet)*fb(2) - sig*bet*fb(3) +sig*fb(1) +fy(2);
    X2_Preston(t) = kapp*alph*bet*fa(2) + (1-alph)*bet*fa(1);
end

eq18_preston = -z_Preston(2,:) -sig*z_Preston(3,:) + X1_Preston;
eq19_preston = -z_Preston(1,:) + kapp*z_Preston(2,:) + X2_Preston;
eq27_preston = -z_Preston(3,:) + psi_pi*z_Preston(1,:) + psi_x*z_Preston(2,:) + s(1,:);

% Peter
X1_peter = zeros(1,T);
X2_peter = zeros(1,T);
fb = zeros(k,T);
fa = zeros(k,T);
fy = zeros(n,T);
for t=1:T
    if t==1
        at_1 = phi0(1,:)';
        bt_1 = phi0(2:end,:)';
        fb(:,t) = at_1/(1-bet) + bt_1*(eye(n)-bet*C)^(-1)*s(:,t);
        fa(:,t) = at_1/(1-alph*bet) + bt_1*(eye(n)-alph*bet*C)^(-1)*s(:,t);
        [faa, fbb, fyy] = fafbfy(at_1,bt_1,param,s(:,t));
    else
        %         fb(:,t) = a_Peter(:,t-1)/(1-bet) + squeeze(b_Peter(:,:,t-1))*(eye(n)-bet*C)^(-1)*s(:,t);
        %         fa(:,t) = a_Peter(:,t-1)/(1-alph*bet) + squeeze(b_Peter(:,:,t-1))*(eye(n)-alph*bet*C)^(-1)*s(:,t);
        fb(:,t) = a_Peter(:,t-1)/(1-bet) + squeeze(b_Peter(:,:,t-1))/(eye(n)-bet*C)*s(:,t);
        fa(:,t) = a_Peter(:,t-1)/(1-alph*bet) + squeeze(b_Peter(:,:,t-1))/(eye(n)-alph*bet*C)*s(:,t);
        [faa, fbb, fyy] = fafbfy(a_Peter(:,t-1),squeeze(b_Peter(:,:,t-1)),param,s(:,t));
    end
    fy(:,t) = f'*(eye(n)-bet*C)^(-1)*s(:,t);
    diff = max(max(abs([fa(:,t)' - faa', fb(:,t)'-fbb', fy(:,t)'-fyy'])));
    if diff > 1e-10
        warning('fs arent the same')
    end
    % these are the LR expectations on RHS of 18 and 19:
        X1_peter(t) = (1-bet)*fb(2,t) - sig*bet*fb(3,t) +sig*fb(1,t) +fy(2,t);
        X2_peter(t) = kapp*alph*bet*fa(2,t) + (1-alph)*bet*fa(1,t);
%     X1_peter(t) = (1-bet)*FB(2,t) - sig*bet*FB(3,t) +sig*FB(1,t) +FY(2,t);
%     X2_peter(t) = kapp*alph*bet*FA(2,t) + (1-alph)*bet*FA(1,t);
end
eq18_peter = -z_Peter(2,:) -sig*z_Peter(3,:) + X1_peter;
eq19_peter = -z_Peter(1,:) + kapp*z_Peter(2,:) + X2_peter;
eq27_peter = -z_Peter(3,:) + psi_pi*z_Peter(1,:) + psi_x*z_Peter(2,:) + s(1,:);

truncated_expectations =0;
if truncated_expectations ==1
    % Construct alternative, truncated long-horizon expectations:
    H=1000; % you need a sum of fcsts up to a min of H=1000 for asymptotics to kick in
    some_t = 4;
    [X1_alt_st, X2_alt_st] = lrexp(H,a_Peter(:,some_t-1),squeeze(b_Peter(:,:,some_t-1)),param,s(:,some_t))
    [X1_alt2_st, X2_alt2_st] = lrexp2(fa(:,some_t),fb(:,some_t),fy(:,some_t),param)
    [X1_peter(some_t), X2_peter(some_t)]
    
    [fa_straight_st, fa_queer_st]=fatest(H,a_Peter(:,some_t-1),squeeze(b_Peter(:,:,some_t-1)),param,s(:,some_t))
    [fb_straight_st, fb_queer_st]=fbtest(H,a_Peter(:,some_t-1),squeeze(b_Peter(:,:,some_t-1)),param,s(:,some_t))
    
    % now construct them for each t
    X1_alt = zeros(1,T);
    X1_alt2 = zeros(1,T);
    X2_alt = zeros(1,T);
    X2_alt2 = zeros(1,T);
    fa_straight = zeros(k,T);
    fa_queer = zeros(k,T);
    fb_straight = zeros(k,T);
    fb_queer = zeros(k,T);
    tic
    for tt=2:T
        [X1_alt(tt), X2_alt(tt)] = lrexp(H,a_Peter(:,tt-1),squeeze(b_Peter(:,:,tt-1)),param,s(:,tt));
        [X1_alt2(tt), X2_alt2(tt)] = lrexp2(fa(:,tt),fb(:,tt),fy(:,tt),param); % uses fa, fb, fy from above explicitly
        [fb_straight(:,tt), fb_queer(:,tt)]=fbtest(H,a_Peter(:,tt-1),squeeze(b_Peter(:,:,tt-1)),param,s(:,tt));
        [fa_straight(:,tt), fa_queer(:,tt)]=fatest(H,a_Peter(:,tt-1),squeeze(b_Peter(:,:,tt-1)),param,s(:,tt));
    end
    toc
    
    until_when =100;
    figure
    set(gcf,'color','w'); % sets white background color
    set(gcf, 'Position', get(0, 'Screensize')); % sets the figure fullscreen
    plot(X1_peter(1:until_when)); hold on
    plot(X1_alt(1:until_when), '*')
    plot(X1_alt2(1:until_when), '--')
    legend('X1 Peter', 'X1 alt 1', 'X1 alt 2')
    
    until_when =100;
    figure
    set(gcf,'color','w'); % sets white background color
    set(gcf, 'Position', get(0, 'Screensize')); % sets the figure fullscreen
    plot(fb_straight(3,1:until_when)); hold on
    plot(fb_queer(3,1:until_when), '--')
    legend('fb straight', 'fb queer')
    
    figure
    set(gcf,'color','w'); % sets white background color
    set(gcf, 'Position', get(0, 'Screensize')); % sets the figure fullscreen
    plot(fa_straight(3,1:until_when)); hold on
    plot(fa_queer(3,1:until_when), '--')
    legend('fa straight', 'fa queer')
    
    eq18_peter_alt = -z_Peter(2,:) -sig*z_Peter(3,:) + X1_alt; % these now equal the previous ones
    eq19_peter_alt = -z_Peter(1,:) + kapp*z_Peter(2,:) + X2_alt;
    figure
    set(gcf,'color','w'); % sets white background color
    set(gcf, 'Position', get(0, 'Screensize')); % sets the figure fullscreen
    plot(eq18_peter_alt); hold on
    plot(eq19_peter_alt, '*')
    
    eq18_zz = -zz(2,:) -sig*zz(3,:) + X1_peter; % also these are the same as before
    eq19_zz = -zz(1,:) + kapp*zz(2,:) + X2_peter;
end



%% Correlations
RHO = zeros(k,k); % 1st dim is per variable, 2nd per model pair
for i=1:k
    % Inflation
    eta = cov(z_BM(i,:), z_Preston(i,:));
    RHO(i,1) = eta(1,2) / (sqrt(eta(1,1))*sqrt(eta(2,2)));
    eta = cov(z_BM(i,:), z_Peter(i,:));
    RHO(i,2) = eta(1,2) / (sqrt(eta(1,1))*sqrt(eta(2,2)));
    eta = cov(z_Peter(i,:), z_Preston(i,:));
    RHO(i,3) = eta(1,2) / (sqrt(eta(1,1))*sqrt(eta(2,2)));
end

corr_table = RHO;
if output_table==1
    columnLabels = {'BM - Preston','BM - Peter', 'Preston - Peter'};
    rowLabels = {'$\pi$', 'Ouput gap', 'Nominal int. rate'};
    matrix2latex_black(corr_table, [tablepath,'/materials11_corrtable.tex'], 'rowLabels', rowLabels, 'columnLabels', columnLabels, ...
        'alignment', 'c', 'format', '%-6.3f', 'size', 'small');
end

%% Plots

learnin=1; % for now, but the idea is to plot the economy only once learning has converged
skipping_step = 1000;
figure
set(gcf,'color','w'); % sets white background color
set(gcf, 'Position', get(0, 'Screensize')); % sets the figure fullscreen
subplot(1,3,1)
plot(z_BM(1,learnin:skipping_step:end), 'k', 'linewidth',2); hold on
plot(z_BM(2,learnin:skipping_step:end), 'b', 'linewidth',2);
plot(z_BM(3,learnin:skipping_step:end), 'r', 'linewidth',2);
ax = gca; % current axes
ax.FontSize = fs;
grid on
grid minor
legend('\pi', 'x', 'i')
xlabel('time')
% ax.YAxis.Exponent = 0;
title('B & M economy')

subplot(1,3,2)
plot(z_Preston(1,learnin:skipping_step:end), 'k', 'linewidth',2); hold on
plot(z_Preston(2,learnin:skipping_step:end), 'b', 'linewidth',2);
plot(z_Preston(3,learnin:skipping_step:end), 'r', 'linewidth',2);
ax = gca; % current axes
ax.FontSize = fs;
grid on
grid minor
legend('\pi', 'x', 'i')
xlabel('time')
% ax.YAxis.Exponent = 0;
title('P economy, P''s ALM')

subplot(1,3,3)
plot(z_Peter(1,learnin:skipping_step:end), 'k', 'linewidth',2); hold on
plot(z_Peter(2,learnin:skipping_step:end), 'b', 'linewidth',2);
plot(z_Peter(3,learnin:skipping_step:end), 'r', 'linewidth',2);
ax = gca; % current axes
ax.FontSize = fs;
grid on
grid minor
legend('\pi', 'x', 'i')
xlabel('time')
% ax.YAxis.Exponent = 0;
title('P economy, Peter''s ALM')
if print_figs ==1
    figname = ['materials11_preston_states']
    cd(figpath)
    export_fig(figname)
    cd(current_dir)
    close
end

% Let's do the same plot, by variables
figure
set(gcf,'color','w'); % sets white background color
set(gcf, 'Position', get(0, 'Screensize')); % sets the figure fullscreen
subplot(1,3,1)
plot(z_BM(1,learnin:skipping_step:end), 'k', 'linewidth',2); hold on
plot(z_Preston(1,learnin:skipping_step:end), 'b', 'linewidth',2);
plot(z_Peter(1,learnin:skipping_step:end), 'r', 'linewidth',2);
ax = gca; % current axes
ax.FontSize = fs;
grid on
grid minor
legend('BM', 'Preston', 'Peter')
xlabel('time')
title('Inflation')

subplot(1,3,2)
plot(z_BM(2,learnin:skipping_step:end), 'k', 'linewidth',2); hold on
plot(z_Preston(2,learnin:skipping_step:end), 'b', 'linewidth',2);
plot(z_Peter(2,learnin:skipping_step:end), 'r', 'linewidth',2);
ax = gca; % current axes
ax.FontSize = fs;
grid on
grid minor
legend('BM', 'Preston', 'Peter')
xlabel('time')
title('Output gap')

subplot(1,3,3)
plot(z_BM(3,learnin:skipping_step:end), 'k', 'linewidth',2); hold on
plot(z_Preston(3,learnin:skipping_step:end), 'b', 'linewidth',2);
plot(z_Peter(3,learnin:skipping_step:end), 'r', 'linewidth',2);
ax = gca; % current axes
ax.FontSize = fs;
grid on
grid minor
legend('BM', 'Preston', 'Peter')
xlabel('time')
% ax.YAxis.Exponent = 0;
title('Nom. int. rate')
if print_figs ==1
    figname = ['materials11_preston_states2']
    cd(figpath)
    export_fig(figname)
    cd(current_dir)
    close
end

% Convergence of learning
step=100;
selected_diff_BM = diff_BM(1:step:end,:);
selected_diff_Preston = diff_Preston(1:step:end,:);
selected_diff_Peter = diff_Peter(1:step:end,:);
selected_time = 1:step:T;
figure
set(gcf,'color','w'); % sets white background color
set(gcf, 'Position', get(0, 'Screensize')); % sets the figure fullscreen
dbm=plot(selected_time,selected_diff_BM, 'k', 'linewidth',2); hold on
dpreston=plot(selected_time,selected_diff_Preston, 'b', 'linewidth',2);
dpeter=plot(selected_time,selected_diff_Peter, 'r', 'linewidth',2);
ax = gca; % current axes
ax.FontSize = fs;
ylim([0, 1e-6])
grid on
grid minor
title('Max absolute difference in regression coefficients')
legend([dbm, dpreston, dpeter], 'BM', 'Preston', 'Peter')
current_dir = pwd;
if print_figs ==1
    figname = ['materials11_convergence_all']
    cd(figpath)
    export_fig(figname)
    cd(current_dir)
    close
end

close all

% Plot errors in simulation
figure
set(gcf,'color','w'); % sets white background color
set(gcf, 'Position', get(0, 'Screensize')); % sets the figure fullscreen
subplot(3,3,1)
plot(eq13)
ax = gca; % current axes
ax.FontSize = fs;
title('13')
subplot(3,3,2)
plot(eq14)
ax = gca; % current axes
ax.FontSize = fs;
title('14')
subplot(3,3,3)
plot(eq27)
ax = gca; % current axes
ax.FontSize = fs;
title('27')

subplot(3,3,4)
plot(eq18_preston)
ax = gca; % current axes
ax.FontSize = fs;
title('18 Preston')
subplot(3,3,5)
plot(eq19_preston)
ax = gca; % current axes
ax.FontSize = fs;
title('19 Preston')
subplot(3,3,6)
plot(eq27_preston)
ax = gca; % current axes
ax.FontSize = fs;
title('27 Preston')

subplot(3,3,7)
plot(eq18_peter)
ax = gca; % current axes
ax.FontSize = fs;
% ylim([-4,4])
title('18 Peter')
subplot(3,3,8)
plot(eq19_peter)
ax = gca; % current axes
ax.FontSize = fs;
title('19 Peter')
subplot(3,3,9)
plot(eq27_peter)
ax = gca; % current axes
ax.FontSize = fs;
title('27 Peter')
