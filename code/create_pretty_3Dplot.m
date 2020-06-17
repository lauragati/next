function create_pretty_3Dplot(z,xxgrid,yygrid,xlab,ylab,zlab,figname,print_figs)
% create a single 3D plot with of 2 variables (no subplots here!)

% Paths
[current_dir, ~, ~,~,~,figpath] = add_paths;
% Plot configs
[fs] = plot_configs;

ng = length(xxgrid);

figure
set(gcf,'color','w'); % sets white background color
set(gcf, 'Position', get(0, 'Screensize')); % sets the figure fullscreen

mesh(xxgrid, yygrid,reshape(z,[ng,ng]))
xlabel(xlab,'interpreter', 'latex', 'fontsize', fs)
ylabel(ylab,'interpreter', 'latex', 'fontsize', fs)
zlabel(zlab,'interpreter', 'latex', 'fontsize', fs, 'rotation',0)

ax = gca; % current axes
ax.FontSize = fs;
set(gca,'TickLabelInterpreter', 'latex');
grid on
grid minor

if print_figs ==1
    disp(figname)
    cd(figpath)
    export_fig(figname)
    cd(current_dir)
    close
end