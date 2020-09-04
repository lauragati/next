function create_pretty_plot_x_holdon(x,y,legendentries,xlab, ylab,xlmult, ylmult, figname,print_figs)
% create a single plot with several series on it (no subplots here!)
% y = ny x T
[ny,T] = size(y);

this_code = mfilename;
max_no_inputs = nargin(this_code);
if nargin < max_no_inputs
    print_figs = 0;
end

% Paths
[current_dir, ~, ~,~,~,figpath] = add_paths;
% Plot configs
[fs, lw] = plot_configs;

figure
set(gcf,'color','w'); % sets white background color
set(gcf, 'Position', get(0, 'Screensize')); % sets the figure fullscreen
for i=1:ny
    h(i) = plot(x,y(i,:), 'linewidth', lw); hold on
end
% plot(zeros(1,T), 'k--', 'linewidth',lw)
ax = gca; % current axes
ax.FontSize = fs;
set(gca,'TickLabelInterpreter', 'latex');
grid on
grid minor
legend(h, legendentries, 'location', 'southoutside', 'interpreter', 'latex')
legend('boxoff')
ax.XAxis.Limits = round(x([1,end]),1); % this is like the 'keeplimits' for dateticks
if nargin == max_no_inputs
    xl = xlabel(xlab,'interpreter', 'latex');
    xl.Position(1) = xlmult* abs(xl.Position(1));
    
    yl = ylabel(ylab,'interpreter', 'latex');
    yl.Rotation = 0; % rotate
    yl.Position(1) = ylmult(1) * yl.Position(1); % move left
    yl.Position(2) = ylmult(2) * yl.Position(2); % move up
end

if print_figs ==1
    disp(figname)
    cd(figpath)
    export_fig(figname)
    cd(current_dir)
    close
end