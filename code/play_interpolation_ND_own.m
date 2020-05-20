% play_interpolation_ND_own
% A playground for N-dimensional interpolation using own codes and Ryan's
% ndim_simplex.m

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

%% The leading example with csapi - the "truth"

true_fun  = @(x1,x2) x1+x2;
a=0;
b=1;
c=2;
d=3;

% Evaluation points and true f at evaluation points
T=10;
Xe = linspace(a,b,T);
Ye = linspace(c,d,T);
% using the notation that XX is the matrix version of the vector grid
[XXe,YYe] = ndgrid(Xe,Ye);
f = true_fun(XXe,YYe);

% Interpolation points
N = 3; % number of nodes
n = T/(N-1);
Xq = [Xe(1:n:end), Xe(end)];
Yq = [Ye(1:n:end), Ye(end)];

[XXq, YYq]=ndgrid(Xq,Yq);
% evaluate function at query points
fq = true_fun(XXq,YYq);
% evaluate spline at query points
pp = csapi({Xq,Yq},fq);
fhat3 = fnval(pp,{Xe,Ye});

figure
subplot(1,2,1)
mesh(XXe,YYe,f)
title('Approximand')
subplot(1,2,2)
mesh(Xe,Ye,fhat3)
title('Spline with "csapi"')
sgtitle('Original approximand with Matlab-internal approximation')

%% ndim_simplex it - not yet workin
x =cell(1,1);
x{1} = Xq;
x{2} = Yq;
xx = [Xe;Ye]; % Ndim x T evaluation nodes
[alph,Q,fracin] = ndim_simplex(x,xx,f)


%% Chebyshev on 2D
close all

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
[XXq,YYq] = ndgrid(x,y);
w = true_fun(XXq,YYq); 

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Step 4. Compute Chebyshev coefficients
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% might not have to compute them painfully
T = chebyshev(z,n);
% ai = T\w/T; % ??????????
ai = (T*T)\w;

% Then the interpolant is
% 1. Convert the evaluation points to [-1,1] nodes
Z1 =2*((Xe-a)./(b-a)-1);
Z2 =2*((Ye-c)./(d-c)-1);
% 2. Evaluate the Chebys there
T1 = chebyshev(Z1,n);
T2 = chebyshev(Z2,n);
% 3. Interpolant
% pxy = T1*ai*T2'; % ??????????
pxy = T2*(T1*ai)';


figure
subplot(1,2,1)
mesh(XXe,YYe,f)
subplot(1,2,2)
mesh(Xe,Ye,pxy)
sgtitle('2D Chebyshev approximation - not working, as you can see')



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