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


fig = figure;
set(gcf,'color','w'); % sets white background color
set(gcf, 'Position', get(0, 'Screensize')); % sets the figure fullscreen
for i=1:ny
    subplot(1,ny,i)
    plot(zeros(1,T), 'k', 'linewidth',lw); hold on
    h1(i) =   plot(y1(i,:), 'linewidth', lw);  
    h2(i) =   plot(y2(i,:), 'linewidth', lw, 'linestyle', '--');
    
    ax = gca; % current axes
    ax.FontSize = fs/1.2;
    set(gca,'TickLabelInterpreter', 'latex');
    grid on
    grid minor
    title(seriesnames{i}, 'interpreter', 'latex', 'fontsize', fs)
%     legend([h1(i), h2(i)], legendentries, 'location', 'southoutside', 'interpreter', 'latex')
%     legend('boxoff')
    if nargin == max_no_inputs
    xl = xlabel(xlab,'interpreter', 'latex', 'fontsize', fs/2);
    xl.Position(1) = xplus + xl.Position(1);
    end
%     ax.YAxis.Exponent = 0;
%     ax.YRuler.Exponent = 0; % turns off scientific notation
end

% Overall legend
% add a bit space to the figure
% fig = gcf;
% fig.Position(2) = fig.Position(2)+250;
fig.InnerPosition = 0.9*fig.OuterPosition;
% add legend
Lgnd = legend('show',[h1(i), h2(i)], legendentries, 'location', 'southoutside', 'interpreter', 'latex', 'NumColumns', 2);
legend('boxoff')
Lgnd.Position(1) = 0.43;
Lgnd.Position(2) = -0;


if print_figs ==1
    disp(figname)
    cd(figpath)
    export_fig(figname)
    cd(current_dir)
    close
end