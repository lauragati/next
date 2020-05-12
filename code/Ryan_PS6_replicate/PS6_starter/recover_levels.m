function [longsim,longsim_lev] = recover_levels(longsim,param,set,modl)

logidx = [modl.ylog,length(modl.Y)+modl.xlog];
ss = model_ss(param,set);
nt = size(longsim,2);

%Add in SS
longsim(logidx,:) = longsim(logidx,:) + repmat(log(ss(logidx))',[1,nt]);
longsim(logidx==0,:) =longsim(logidx==0,:) + repmat(ss(logidx==0)',[1,nt]);

%Put log-vars into levels
longsim_lev = longsim;
longsim_lev(logidx,:) = exp(longsim_lev(logidx,:));