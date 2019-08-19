function [fk, fpibar, Apibar, fxi, Axi, Sxi, SIG, h0, hpibar, H, R, l] = cemp_functions_matrices(param, k_t_1, pibar_t_1, observation_block_size)
% 22 June 2019
% observation_block_size allows you to choose the size of the obseration (1=the sample with lots of missing values, 2 = the sample with some missings, 3 = sample without missings)
% block based on the availability of data
%
% All of this given particle i

[fk, fpibar, Apibar, fxi, Axi] = cemp_state_block(param, k_t_1, pibar_t_1);
[h0, hpibar, H, R, l] = cemp_obs_block(param,observation_block_size);
[SIG, Sxi] = cemp_SIG_S(param); 

