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
%     ax.Position(1) = ax.Position(1)+0.1; % bottom left x
    ax.Position(2) = ax.Position(2)+0.1; % bottom left y
    ax.Position(3) = 0.9*ax.Position(3); % width
    ax.Position(4) = 0.8*ax.Position(4); % height
    ax.FontSize = fs/1.5;
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
    ax.YAxis.Exponent = 0;
    ax.YRuler.Exponent = 0; % turns off scientific notation
end

% Overall legend
Lgnd = legend('show',[h1(i), h2(i)], legendentries, 'fontsize',fs*0.75,'location', 'southoutside', 'interpreter', 'latex', 'NumColumns', 2);
legend('boxoff')
Lgnd.Position(1) = 0.4; % IRFs
Lgnd.Position(2) = 0.05;



if print_figs ==1
    disp(figname)
    cd(figpath)
    export_fig(figname)
    cd(current_dir)
    close
end