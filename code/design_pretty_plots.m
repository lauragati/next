% design_pretty_plots.m
% Experiment with a unified visual identity for all plots in drafts and
% presentations
print_figs = 0;

% Let's illustrate the visual identity first:
series1 = rand(1,10);
series2 = rand(1,10);
series3 = rand(1,10);

series1b = rand(1,10);
series2b = rand(1,10);
series3b = rand(1,10);

% Note that you have to toggle a zero line manually within the
% plotting-functions.
figname = 'visual_identity1';
create_pretty_plot_x(1:10, series1, figname, print_figs)
figname = 'visual_identity2';
create_pretty_plot_holdon([series1;series2;series3],{'series1', 'series2', 'series3'}, figname, print_figs)
figname = 'visual_identity3';
create_pretty_subplots([series1;series2;series3],{'series1', 'series2', 'series3'}, figname, print_figs)
figname = 'visual_identity4';
create_pretty_subplots_holdon([series1;series2;series3],[series1b;series2b;series3b],{'series1', 'series2', 'series3'},{'Anchoring', 'RE'}, figname, print_figs)


% Figure 1. Market-based inflation expectations -->
% create_motivation_plots.m

% Figure 2. Comparative statics. Policy function. --> analyze_opt_policy.m

% Figure 3. Taylor rule versus optimal policy in the simulation -->
% command_pea.m for Taylor rule and compare_value_pea_results.m for PEA sol

% Figure 4. Central bank loss function as a function of psi_pi -->
% materials16.m, plot_sim_loss.m

% Figures 5,6,7 IRFs after contractionary monpol shock, anchored and unanchored
% --> command_IRFs_anchoring_pretty.m

% Figure 8. Policy function for PEA vs. VFI -->
% compare_value_pea_results.m
