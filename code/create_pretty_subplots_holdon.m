function create_pretty_subplots_holdon(y1,y2,seriesnames,legendentries,xlab,xplus,figname,print_figs)
% y = ny x T
[ny,T] = size(y1);

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
    subplot(1,ny,i)
    h1(i) =   plot(y1(i,:), 'linewidth', lw); hold on
    h2(i) =   plot(y2(i,:), 'linewidth', lw);
    plot(zeros(1,T), 'k--', 'linewidth',lw)
    ax = gca; % current axes
    ax.FontSize = fs/2;
    set(gca,'TickLabelInterpreter', 'latex');
    grid on
    grid minor
    title(seriesnames{i}, 'interpreter', 'latex', 'fontsize', fs)
    legend([h1(i), h2(i)], legendentries, 'location', 'southoutside', 'interpreter', 'latex')
    legend('boxoff')
    if nargin == max_no_inputs
    xl = xlabel(xlab,'interpreter', 'latex', 'fontsize', fs/2.5);
    xl.Position(1) = xplus + xl.Position(1);
    end
end


if print_figs ==1
    disp(figname)
    cd(figpath)
    export_fig(figname)
    cd(current_dir)
    close
end