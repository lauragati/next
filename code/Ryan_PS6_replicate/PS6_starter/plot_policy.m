function f = plot_policy(xgrid,pol_init,pol_final)

nk = length(xgrid{1});
na = length(xgrid{2});

kpts = [1 ceil(nk/2), nk];
apts = [1 ceil(na/2), na];


%Plot initial and final policy functions together
f = figure;

for jj = 1:3
    %With respect to S, for two different news shocks
    subplot(1,2,1);hold on;

    xxgrid = [xgrid{1};repmat(xgrid{2}(apts(jj)),[1,nk])];
    hh1 = ndim_simplex_eval(xgrid,xxgrid,pol_init);
    hh2 = ndim_simplex_eval(xgrid,xxgrid,pol_final);
    
    plot(xgrid{1},hh1, '-bo');
    plot(xgrid{1},hh2, '-gs');
    
    xlabel('Capital');
    
    subplot(1,2,2);hold on;
    xxgrid = [repmat_col(xgrid{1}(kpts(jj)),na);xgrid{2}];
    hh3 = ndim_simplex_eval(xgrid,xxgrid,pol_init);
    hh4 = ndim_simplex_eval(xgrid,xxgrid,pol_final);
    
    plot(xgrid{2},hh3, '-bo');
    plot(xgrid{2},hh4, '-gs');
    
    
    xlabel('Productivity');
    
end


%Title and legend
legend('Linearized Model','Nonlinear (Collocation)')

