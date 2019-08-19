function [SIG, Sxi] = cemp_SIG_S(param) 
sige  = param(7);
sigmu = param(8);

Sxi = [1,1; 1,0; 1,1]; % shock impact in lin state equation % <--
% Sxi = [1,1; 1,1; 1,1]; 

SIG = diag([sige sigmu]); %  VC of state innovations.