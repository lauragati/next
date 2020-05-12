% PARAMETETRS - This function returns a parameter structure to use in the model solution.
%
%


function [param,set] = parameters()



%Parameters are values that are to be estiamted. They are passed as the 
%first argument to model_prog.m.
param.sige = .1; %Standard deviation of TFP growths shocks


%Settings are values that you 'calibrate' and therefore do change in
%estimation. They are passed as the second argument to model_prog.m
set.bet = .98;  %Never call anything beta...it's a matlab function
set.chi = 1;    %Disutility of labor
set.del = .05;  %Depreciation Rate
set.alph = .36; %Capital Share (never call anything alpha, either)
set.phi = 1;    %Irreversible investment parameter
set.gam = 1.01;
set.kbar = NaN;

%Choose whether to consider measurement error or a higher-order 
%approximation.
set.me_eq = [];      %Empty vector means no measurement error
set.approx_deg = 1;  %Degree of Approximation