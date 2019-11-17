function FEV = obj_minFEV(param,gx,hx,SIG,T,burnin,e,Aa,Ab,As)

[~, paramset] = parameters_next;
fn = fieldnames(paramset);
make_index(fn);

% create one struct with all the stuff for sim_learn and its children
paramset.(fn{gbar_idx}) = param; % in this struct, replace gbar

PLM      = 1; % constant only
cgain    = 3; % constant gain
critCEMP = 1; % this doesn't matter

[~, ~, ~, ~, ~, ~, FEt_1_c] = sim_learn(gx,hx,SIG,T,burnin,e, Aa, Ab, As, paramset, PLM, cgain, critCEMP);

FEV = nanmean(FEt_1_c.^2); %FEV is integrating over time for a given history
