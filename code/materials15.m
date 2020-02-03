% materials15

% Goals:
% 1.) understand difference in  behavior CEMP vs. CUSUM criterion
% 29 Jan 2020
clearvars
close all
clc

% Add all the relevant paths and grab the codename
this_code = mfilename;
[current_dir, basepath, BC_researchpath,toolpath,export_figpath,figpath,tablepath,datapath] = add_paths;

% Variable stuff ---
print_figs        = 0;
stop_before_plots = 0;
skip_old_plots    = 0;
output_table = print_figs;

skip_old_stuff = 1;

%% Do IRFs, simulation and gains for anchoring model, CEMP or CUSUM criterion
command_IRFs_anchoring

%% Do fmincon for psi_pi in the anchoring model
grid_search

%% Compare phi_pi* and psi_pi* for the case when all shocks have a common rho
rhoall = 0.3;
lami=1;
[psi_pi_opt, phi_pi_opt] =  optTR_shared_rho(rhoall,lami);

lami_values = linspace(0,1);
rho_values  = linspace(0,1);
M =size(rho_values,2);
psi = zeros(M,M);
phi = zeros(M,M);
for i=1:M
    lami=lami_values(i);
    for j=1:M
        rhoall=rho_values(j);
        [psi(i,j), phi(i,j)] = optTR_shared_rho(rhoall,lami);
    end
end

%% Plots
yaxislimits= [-0.1, 5];
% Fix lami
figure
plot(rho_values, psi(2,:)); hold on
plot(rho_values, phi(2,:), '--')
legend('Learning', 'RE')
xlabel('\rho')
ylim(yaxislimits)
xlim([rho_values(1), rho_values(end)])
sgtitle(['\lambda_i = ', num2str(lami_values(2))])

figure
plot(rho_values, psi(50,:)); hold on
plot(rho_values, phi(50,:), '--')
legend('Learning', 'RE')
xlabel('\rho')
ylim(yaxislimits)
xlim([rho_values(1), rho_values(end)])
sgtitle(['\lambda_i = ', num2str(lami_values(50))])


figure
plot(rho_values, psi(end-1,:)); hold on
plot(rho_values, phi(end-1,:), '--')
legend('Learning', 'RE')
xlabel('\rho')
ylim(yaxislimits)
xlim([rho_values(1), rho_values(end)])
sgtitle(['\lambda_i = ', num2str(lami_values(end-1))])


%% Fix rho
figure
plot(lami_values, psi(:,2)); hold on
plot(lami_values, phi(:,2), '--')
legend('Learning', 'RE')
xlabel('\lambda_i')
% ylim(yaxislimits)
xlim([rho_values(1), rho_values(end)])
sgtitle(['\rho = ', num2str(rho_values(2))])

figure
plot(lami_values, psi(:,50)); hold on
plot(lami_values, phi(:,50), '--')
legend('Learning', 'RE')
xlabel('\lambda_i')
% ylim(yaxislimits)
xlim([rho_values(1), rho_values(end)])
sgtitle(['\rho = ', num2str(rho_values(50))])

figure
plot(lami_values, psi(:,end-1)); hold on
plot(lami_values, phi(:,end-1), '--')
legend('Learning', 'RE')
xlabel('\lambda_i')
% ylim(yaxislimits)
xlim([rho_values(1), rho_values(end)])
sgtitle(['\rho = ', num2str(rho_values(end-1))])
