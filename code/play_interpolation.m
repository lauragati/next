% play_interpolation
% A playground for interpolation, linear and spline, and shape-preserving

clearvars
close all
clc

% Add all the relevant paths and grab the codename
this_code = mfilename;
[current_dir, basepath, BC_researchpath,toolpath,export_figpath,figpath,tablepath,datapath] = add_paths;
todays_date = strrep(datestr(today), '-','_');

print_figs=0;
skip = 1;

% Plot configs
[fs, lw] = plot_configs;
%% One-dimensional interpolation

true_fun  = @(x) x.^3;
true_fun2 = @(x) abs(x);
a=-1;
b=1;

Xqmat = linspace(a,b,20); % like the true sample
% f = true_fun(X);
f = true_fun2(Xqmat);

T = length(f);

% "interp" command
fhat1 = interp(f,2); % this one just stretches it by a factor of 2


% prepare nodes = query points Xq, and associated function values Yq
N = 5; % number of nodes
n = T/(N-1);
Xq = [Xqmat(1:n:end), Xqmat(end)];  % query points, i.e. the nodes themselves
Yq = [f(1:n:end), f(end)]; % function values at query points

% the "interp1" command
fhat2 = interp1(Xq,Yq,Xqmat); % linear interpolation
fhat3 = interp1(Xq,Yq,Xqmat, 'spline'); % cubic spline interpolation
fhat4 = interp1(Xq,Yq,Xqmat, 'pchip'); % shape-preserving cubic hermite interpolating polynomial

% The "spline" command
% this gives you the coeffs, note that pp is a struct
pp = spline(Xq,Yq);
nodes = pp.breaks;
a = pp.coefs;
n_intervals =pp.pieces;
% this gives you the interpolant
fhat5 = spline(Xq,Yq,Xqmat);
% you can also do it using the coeffs explicity
fhat5_alt = ppval(pp,Xqmat);
fhat5-fhat5_alt; % bingo

% The "pchip" command - works exactly like the spline command
ppp= pchip(Xq,Yq);
fhat4_alt = ppval(ppp,Xqmat);
fhat4-fhat4_alt; % bingo

% Comparing pchip with SPLINE:
%     The function s(x) supplied by SPLINE is constructed in exactly the same way,
%     except that the slopes at the X(j) are chosen differently, namely to make
%     even D^2s(x) continuous. This has the following effects.
%     SPLINE is smoother, i.e., D^2s(x) is continuous.
%     SPLINE is more accurate if the data are values of a smooth function.
%     pchip has no overshoots and less oscillation if the data are not smooth.
%     pchip is less expensive to set up.
%     The two are equally expensive to evaluate.


% let's compare with my thing
% set up objective spline with equation system
resids = obj_spline_natural(a,Xq,Yq);

% %Optimization Parameters
options = optimoptions('fsolve', 'TolFun', 1e-9, 'display', 'iter', 'MaxFunEvals', 10000);
tic
objh = @(coeffs) obj_spline_natural(coeffs,Xq,Yq);
[coeffs_opt,FVAL] = fsolve(objh,a, options);
toc

a-coeffs_opt % there's the problem
% although Matlab's spline uses the "not-a'knot" end conditions (specifies
% third derivates at endnodes) and not a natural spline, I have reason to
% suspect that my spline is not quite doing the right thing.

%% Plot

figure
set(gcf,'color','w'); % sets white background color
set(gcf, 'Position', get(0, 'Screensize')); % sets the figure fullscreen
plot(Xqmat,f, 'linewidth',lw); hold on
plot(Xqmat, fhat2, 'linewidth',lw)
plot(Xqmat, fhat3, '--', 'linewidth',lw)
plot(Xqmat, fhat4, 'o', 'linewidth',lw)
plot(Xqmat, fhat5, 'x', 'linewidth',lw)
legend('Approximand', 'Linear interpolation with "interp1"', 'Cubic spline with "interp1"', 'Shape-preserving cubic with "interp1"', ...
    'Cubic spline with "spline"', ...
    'location', 'northwest')
ax = gca; % current axes
ax.FontSize = fs;
grid on
grid minor

%% FAZIT

% 1.) interp1('spline') = spline
% 2.) The logic in all these is to use a small number of nodes but pass a
% bigger number of data thru 'em
% 3.) Shape-preserving cubic spline thru the pchip option works well when
% there are kinks to avoid spline-like oscillations around 'em
% 4) there is a pchip command too> it's Piecewise Cubic Hermite Interpolating Polynomial


%% A multidimensional problem - a vector-valued interpolant
close all
clearvars
clc

true_fun3  = @(x1,x2) x1.^(1/2)+x2.^(1/2);
a=0;
b=1;

T=100;
X1 = linspace(a,b,T);
X2 = linspace(a,b,T);
[Xqmat,Y] = ndgrid(X1,X2); % I think both meshgrid and ndgrid work fine: it's just the order of the coordinates that's different.
% [X,Y] = meshgrid(X1,X2);

f = true_fun3(Xqmat,Y);


figure
mesh(Xqmat,Y,f)
close

N = 3; % number of nodes
n = T/(N-1);
Xq1 = [X1(1:n:end), X1(end)];
Xq2 = [X2(1:n:end), X2(end)];
[Xq,Yq] = ndgrid(Xq1,Xq2);
V = true_fun3(Xq,Yq);

Vq1 = interpn(Xq,Yq,V,Xqmat,Y);
Vq2 = interpn(Xq,Yq,V,Xqmat,Y,'spline');
% Vq3 = interpn(Xq,Yq,V,X,Y, 'pchip'); % ain't supported in 2D or more

figure
subplot(1,3,1)
mesh(Xqmat,Y,f)
title('Approximand')
subplot(1,3,2)
mesh(Xqmat,Y,Vq1)
title('Linear interpolation')
subplot(1,3,3)
mesh(Xqmat,Y,Vq2)
title('Spline interpolation')

% that's cool, but I'd like to get the coefficients too

%% a scalar-valued interpolant with several input variables
clc

true_fun4  = @(x) x(1,:).^(1/2)+x(2,:).^(1/2);
a=0;
b=1;

T=10;
X1 = linspace(a,b,T);
X2 = linspace(a,b,T);

N = 3; % number of nodes
n = T/(N-1);
Xq1 = [X1(1:n:end), X1(end)];
Xq2 = [X2(1:n:end), X2(end)];

% evaluation point matrix X (the "sample")
[X1,X2] = ndgrid(X1,X2);
[row, col] = size(X1);
Xqmat = zeros(2,row*col);
index=0;
for i=1:row
    for j=1:col
        index=index+1;
        Xqmat(:,index) = [X2(i,j);X1(i,j)]; % put kgrid on top, ggrid below
    end
end

% true function at evaluation points
f = true_fun4(Xqmat);

% query point matrix Xq
[X1,X2] = ndgrid(Xq1,Xq2);
[row, col] = size(X1);
Xq = zeros(2,row*col);
index=0;
for i=1:row
    for j=1:col
        index=index+1;
        Xq(:,index) = [X2(i,j);X1(i,j)]; % put kgrid on top, ggrid below
    end
end

% true function at query points (nodes)
Yq = true_fun4(Xq);

% dodammn darn that doesn't work. X MUST be a vector.
% fhat1 = spline(Xq,Yq,X);

%% Chebyshev approximation in R2
% Judd, Numerical, p. 243 Mac, Algorithm 6.4 (Using the notation there)

true_fun3  = @(x1,x2) x1+x2;
a=0;
b=1;
c=2;
d=3;

T=10;
X1 = linspace(a,b,T);
X2 = linspace(c,d,T);

[X,Y] = ndgrid(X1,X2);
f = true_fun3(X,Y);

m=3; % # of interpolation nodes = order of cheby polyms + 1
n= m-1; % order of cheby polynomials
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Step 1. Compute Chebyshev interpolation nodes on [-1,1]
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
z = -cos((2*(1:m) -1)/(2*m)*pi);

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Step 2. Adjust nodes to [a,b] and [c,d] intervals
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
x = (z+1).*(b-a)/2 +a;
y = (z+1).*(d-c)/2 +c;

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Step 3. Evaluate f at the approximation nodes
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
[Xq1,Xq2] = ndgrid(x,y);
%  [row, col] = size(Xq1);
%  Xqmat = zeros(2,row*col);
%  w =zeros(m,m);
%  index=0;
%  for i=1:row
%      for j=1:col
%          index=index+1;
%          Xqmat(:,index) = [Xq1(i,j);Xq2(i,j)];
%          w_alt(i,j) = true_fun3(x(i),y(j));
%      end
%  end
w = true_fun3(Xq1,Xq2); % --> so you can actually do this really quickly!
% w-w_alt % they are the same

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Step 4. Compute Chebyshev coefficients
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% might not have to compute them painfully
T = chebyshev(z,n);
ai = T\w/T; 

% Then the interpolant is
% 1. Convert the evaluation points to [-1,1] nodes
Z1 =2*((X1-a)./(b-a)-1);
Z2 =2*((X2-c)./(d-c)-1);
% 2. Evaluate the Chebys there
T1 = chebyshev(Z1,n);
T2 = chebyshev(Z2,n);
% 3. Interpolant
pxy = T1*ai*T2';

figure
subplot(1,2,1)
mesh(X,Y,f)
subplot(1,2,2)
mesh(X1,X2,pxy)


function Tx = chebyshev(x,n) % checked - correct
[xm,xn]=size(x);
tau=max(xm,xn); % the number of data points
x=reshape(x,tau,1);
Tx = zeros(tau,n+1);
% creates Chebyshev polynomials for all x, of order 0,...,n
Tx(:,1) = ones(tau,1); % order 0
Tx(:,2) = x.*ones(tau,1);
for j=3:n+1
    Tx(:,j) = 2*x.*Tx(:,j-1)-Tx(:,j-2);
end
end




