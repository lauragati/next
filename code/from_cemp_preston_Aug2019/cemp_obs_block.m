function [h0, hpibar, H, R, l] = cemp_obs_block(param,observation_block_size)
% 26 June 2019
% observation_block_size allows you to choose the size of the obseration (1=the sample with lots of missing values, 2 = the sample with some missings, 3 = sample without missings)
% block based on the availability of data
%
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

% Observation block
Hpirow  = [0,0,0,1];
Hpi1row = [1-gam, 0, rho, gam];
Hpi2row = [1-gam^2, 0 rho*(gam+rho), gam^2];

switch observation_block_size
    case 1
        % 'lots missing'
        bigH = vertcat(Hpirow, Hpi1row);
        l = size(bigH,1);
        hpibar = bigH(:,1);
        H = bigH(:,2:end)'; % note: I'm transposing H so that I can be consistent with CEMP's notation
        
        h0 = pistar*ones(l,1);
        R = diag([sigo1,sigo4]);
        
    case 2
        % 'some missing'
        bigH = vertcat(Hpirow, Hpi1row, Hpi2row, Hpi1row);
        l = size(bigH,1);
        hpibar = bigH(:,1);
        H = bigH(:,2:end)';
        
        h0 = pistar*ones(l,1);
        R = diag([sigo1,sigo2,sigo3,sigo4]);
        
    case 3
        % 'full'
        bigH = vertcat(Hpirow, Hpi1row, Hpi2row, Hpi1row, Hpi1row);
        l = size(bigH,1);
        hpibar = bigH(:,1);
        H = bigH(:,2:end)';
        
        h0 = pistar*ones(l,1);
        R = diag([sigo1,sigo2,sigo3,sigo4,sigo5]);
end