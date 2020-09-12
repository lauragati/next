function create_pretty_plot_x(x,y,xlab,ylab,xlplus, ylplus, figname,print_figs)
% create a single plot with several series on it (no subplots here!)
% y = ny x T
[ny,T] = size(y);

this_code = mfilename;
max_no_inputs = nargin(this_code);
if nargin < max_no_inputs %no shock specified
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
    plot(x,y(i,:), 'linewidth', lw); hold on
end
% plot(x,zeros(1,T), 'k--', 'linewidth',lw)
ax = gca; % current axes
ax.FontSize = fs;
set(gca,'TickLabelInterpreter', 'latex');
grid on
grid minor
if nargin == max_no_inputs
    xl = xlabel(xlab,'interpreter', 'latex', 'fontsize', fs*4/5);
    xl.Position(1) = xlplus(1) + xl.Position(1);
    xl.Position(2) = xlplus(2) + xl.Position(2);

    
    yl = ylabel(ylab,'interpreter', 'latex', 'fontsize', fs*4/5);
    yl.Rotation = 0; % rotate
    yl.Position(1) = ylplus(1) + yl.Position(1); % move left
    yl.Position(2) = ylplus(2) + yl.Position(2); % move up

end


if print_figs ==1
    disp(figname)
    cd(figpath)
    export_fig(figname)
    cd(current_dir)
    close
end