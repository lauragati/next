function h = plot_fe_histograms(FE)
% specific histogram-plotting for forecast errors in estimation
% 20 July 2020

[fs, lw] = plot_configs;

fesim = squeeze(FE(1,:,:));
femean = mean(fesim,2);
fesimmax = max(fesim,[],2);
fesimmin = min(fesim,[],2);

% Plot histograms of fe
figure
set(gcf,'color','w'); % sets white background color
set(gcf, 'Position', get(0, 'Screensize')); % sets the figure fullscreen
subplot(2,2,1)
hist(femean)
ax = gca; % current axes
ax.FontSize = fs*3/4;
set(gca,'TickLabelInterpreter', 'latex');
grid on
grid minor
title('Mean fe in cross-section','interpreter', 'latex', 'fontsize', fs*3/4)
subplot(2,2,2)
hist(fesimmax)
ax = gca; % current axes
ax.FontSize = fs*3/4;
set(gca,'TickLabelInterpreter', 'latex');
grid on
grid minor
title('Max fe across histories','interpreter', 'latex', 'fontsize', fs*3/4)
subplot(2,2,3)
hist(fesimmin)
ax = gca; % current axes
ax.FontSize = fs*3/4;
set(gca,'TickLabelInterpreter', 'latex');
grid on
grid minor
title('Min fe across histories','interpreter', 'latex', 'fontsize', fs*3/4)
subplot(2,2,4)
hist(fesim(:))
ax = gca; % current axes
ax.FontSize = fs*3/4;
set(gca,'TickLabelInterpreter', 'latex');
grid on
grid minor
title('Fe across all cross-sections','interpreter', 'latex', 'fontsize', fs*3/4)
