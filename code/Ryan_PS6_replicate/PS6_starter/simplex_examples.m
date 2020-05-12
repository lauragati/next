% SIMPLEX_EXAMPLES - Code from class showing how to use the linean
% simplex/interpolation code suggested for the problem set.


%EX1: One dimensional example

Xpoints  = 1:10;
Xvals    = Xpoints;    %Initial evaluation points same as grid points
Yvals  = randn(1,10);

%Create the function weights a(i)
Xgrid = cell(1,1);
Xgrid{1} = Xpoints;
alph = ndim_simplex(Xgrid,Xvals,Yvals);

%Now evaluate points off the grid
Xextend = 1:.1:10;
Yoffgrid = ndim_simplex_eval(Xgrid,Xextend,alph);

figure
plot(Xvals,Yvals, '-x','linewidth', 2);hold on;
plot(Xextend,Yoffgrid,'-o');
pause
close

%EX2 - Overdetermination Case
Xpoints = 1:10;
Xvals = 1:.5:10;
Yvals = randn(1,19);

%Create the function weights a(i)
Xgrid = cell(1,1);
Xgrid{1} = Xpoints;
alph = ndim_simplex(Xgrid,Xvals,Yvals);

%Evaluate points ON the grid
Yongrid = ndim_simplex_eval(Xgrid,Xpoints,alph);
figure
plot(Xvals,Yvals, 'x','linewidth', 2);hold on;
plot(Xpoints,Yongrid,'-o');
pause
close

%EX3 - Two dimensions
f = @(x,z) randn(1,length(x));% x.^3 + z.^.25;
Xpoints = 1:.2:2;
Zpoints = 1:2:9;
[XXpoints,ZZpoints] = ndgrid(Xpoints,Zpoints);
XXvals = [XXpoints(:)'; ZZpoints(:)'];
Yvals = f(XXvals(1,:),XXvals(2,:));

%Create the function weights a(i)
Xgrid = cell(1,2);
Xgrid{1} = Xpoints;
Xgrid{2} = Zpoints;
alph = ndim_simplex(Xgrid,XXvals,Yvals);

%Evaluate and plot a more refined
Xex = 1:.01:2;
Zex = 1:.01:9;
[XXex,ZZex] = ndgrid(Xex,Zex);
XXgrid_ex = [XXex(:)'; ZZex(:)'];
fvals = ndim_simplex_eval(Xgrid,XXgrid_ex,alph);

%Plot
figure
mesh(XXex,ZZex,reshape(fvals, size(XXex)))
pause
close

