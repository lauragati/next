%
% mom_tab - make the basic moment tables to screen using simulated data
%
% usage:
% 
% mom_tab(simdat, var_idx, var_names)
%
% where
%
% simdat  = nvar-by-nt matrix of simulated data.
% var_idx = indexies of the relevant variables in Y (first value should index output)
% var_names = names of the variables

function [sig ar_rho rho, SigY0] = mom_tab_sim(simdat, var_idx, var_names, varargin)

nt = length(simdat);

if nargin>3
    prefs = varargin{1};
else
    prefs.ndec = 8;
end

sformat = ['%1.' num2str(prefs.ndec), 'f\t'];


simdat = demean(simdat')';


SigY0 = 1/nt*(simdat*simdat');
SigY1 = 1/(nt-1)*simdat(:,2:end)*simdat(:,1:end-1)';


yidx = var_idx(1);
tit_str = [];
for j = 1:length(var_idx)
    idx = var_idx(j);
    sig(j) = sqrt(SigY0(idx,idx));
    rho(j) = SigY0(idx,yidx)/sqrt(SigY0(idx,idx)*SigY0(yidx,yidx));
    ar_rho(j) = SigY1(idx,idx)/(sig(j)^2);
    tit_str = [tit_str, var_names{j}, '\t'];
end


    
disp('Standard Deviations')
disp(sprintf(tit_str))
disp(sprintf(sformat, 100*sig))
disp(' ');


disp('Stddev(X)/Stddev(Y)')
disp(sprintf(tit_str))
disp(sprintf(sformat, sig/(sig(1))))
disp(' ');


disp('Auto-correlations')
disp(sprintf(tit_str))
disp(sprintf(sformat, ar_rho))
disp(' ');

disp('Correlation w/ Y')
disp(sprintf(tit_str))
disp(sprintf(sformat, rho))
disp(' ');
