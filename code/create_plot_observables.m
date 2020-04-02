function create_plot_observables(y,seriesnames,figtitle,figname,print_figs)
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
        subplot(1,ny,i)
        plot(y(i,:), 'linewidth', lw); hold on
        plot(zeros(1,T), 'k--', 'linewidth',lw)
        ax = gca; % current axes
        ax.FontSize = fs;
        grid on
        grid minor
    end
elseif nargin > 1
    figure
    set(gcf,'color','w'); % sets white background color
    set(gcf, 'Position', get(0, 'Screensize')); % sets the figure fullscreen
    for i=1:ny
        subplot(1,ny,i)
        h(i) =   plot(y(i,:), 'linewidth', lw); hold on
        plot(zeros(1,T), 'k--', 'linewidth',lw)
        ax = gca; % current axes
        ax.FontSize = fs;
        grid on
        grid minor
        title(seriesnames{i})
    end
    sgtitle(figtitle, 'FontSize',fs)
end


if print_figs ==1
    disp(figname)
    cd(figpath)
    export_fig(figname)
    cd(current_dir)
    close
end