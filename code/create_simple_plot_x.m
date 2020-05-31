function create_simple_plot_x(x,y,seriesnames,figtitle,figname,print_figs)
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
linestyles = {'-k','-','--', ':'};

if nargin ==1
    figure
    set(gcf,'color','w'); % sets white background color
    set(gcf, 'Position', get(0, 'Screensize')); % sets the figure fullscreen
    for i=1:ny
        plot(x,y(i,:), 'linewidth', lw); hold on
    end
    plot(x,zeros(1,T), 'k--', 'linewidth',lw)
    ax = gca; % current axes
    ax.FontSize = fs;
    grid on
    grid minor
elseif nargin > 1
    figure
    set(gcf,'color','w'); % sets white background color
    set(gcf, 'Position', get(0, 'Screensize')); % sets the figure fullscreen
    for i=1:ny
        h(i) =  plot(x,y(i,:), 'linewidth', lw); hold on
    end
    plot(x,zeros(1,T), 'k--', 'linewidth',lw)
    ax = gca; % current axes
    ax.FontSize = fs;
    grid on
    grid minor
    legend(h,seriesnames, 'interpreter', 'latex' )
    title(figtitle, 'FontSize',fs)
end


if print_figs ==1
    disp(figname)
    cd(figpath)
    export_fig(figname)
    cd(current_dir)
    close
end