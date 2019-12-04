% empirical_cgain_values.m
% just summarizes values for the constant gain kappa found in the literature:


% estimated
estim.eusepi_giannoni_preston2019limits = 0.05;
estim.CEMP = 0.145;
estim.milani2007 = 0.0183;
estim.branch_evans2006 = 0.062;
% estim.malmendier_nagel2016 = 3.044; % but their def may not be compatible
% with the others


% calibrated
calib.eusepi_preston2011 = 0.002; 
calib.williams2003_1 = 0.1;
calib.williams2003_2 = 0.05;
calib.williams2003_3 = 0.03;
calib.orphanides_williams2005 = 0.02;
% calib.gerko2019_1 = estim.malmendier_nagel2016; This is for symmetric
% gains. It's strange though because it's an order of magnitude larger than
% the numbers she otherwise takes.
calib.gerko2019_2 = 0.01; % but Gerko considers "asymmetric gains:" one gain on inflation, one on the output gap.
calib.gerko2019_3 = 0.02;
calib.gerko2019_1 = 0.05;

% avg estimate
sum_calib=0;
fnc = fieldnames(calib);
for k=1:numel(fnc)
    sum_calib = sum_calib + calib.(fnc{k});
end
avg_calib = sum_calib / numel(fnc);

% avg estimate
sum_estim=0;
fne = fieldnames(estim);
for k=1:numel(fne)
    sum_estim = sum_estim + estim.(fne{k});
end
sum_wo_CEMP = sum_estim - estim.CEMP;

avg_estim = sum_estim / numel(fne);
avg_estim_wo_CEMP = sum_wo_CEMP / (numel(fne)-1);


sum_overall = sum_calib + sum_estim;
sum_wo_CEMP = sum_overall - estim.CEMP;

avg_overall = sum_overall / (numel(fnc) + numel(fne));
avg_overall_wo_CEMP = sum_wo_CEMP / (numel(fnc) + numel(fne)-1);

disp(['Average = ', num2str(avg_overall), '; average w/o CEMP = ', num2str(avg_overall_wo_CEMP), ...
    '; average estimate = ', num2str(avg_estim), '; average estimate w/o CEMP = ', num2str(avg_estim_wo_CEMP)])

% Observations T quarters old receives the weight: (1-kapp)^T (Eusepi & Preston 2011)


% The number of quarters used to form expectations is given by 1/gain. (Milani 2007)
no_quarters_CEMP = 1/estim.CEMP;
no_quarters_woCEMP = 1/avg_overall_wo_CEMP;
no_quarters_max = 1/calib.eusepi_preston2011;