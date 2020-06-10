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
ndim=2;
true_fun  = @(x1,x2) x1+x2;
a=0;
b=1;
c=2;
d=3;

% Evaluation points and true f at evaluation points
Tx = 120; % # of X-values
Ty = 120; % # of Y-values
T = Tx*Ty; % sample size is dim1*dim2
Xe = linspace(a,b,Tx);
Ye = linspace(c,d,Ty);
% using the notation that XX is the matrix version of the vector grid
[XXe,YYe] = ndgrid(Xe,Ye);
f = true_fun(XXe,YYe);

% Interpolation points
Nx = 25; % number of nodes in the x-dimension
Ny = 5; % number of nodes in the y-dimension
Nq = Nx*Ny; % length of query sample
nx = Tx/(Nx-1);
ny = Ty/(Ny-1);
Xq = [Xe(1:nx:end), Xe(end)];
Yq = [Ye(1:ny:end), Ye(end)];

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
mesh(XXe,YYe,fhat3)
title('Spline with "csapi"')
sgtitle('Original approximand with Matlab-internal approximation')

%% Let's try spapi, my hope for shape-preserving ND interpolation
% For spapi, you can define the order of the interpolating polynomial in
% each dimension (first cell argument)
pp = spapi({3,3},{Xq,Yq},fq);
fhat4 = fnval(pp,{Xe,Ye});

% let's do a piece-wise linear
pp = spapi({2,2},{Xq,Yq},fq);
fhat5 = fnval(pp,{Xe,Ye});

% create spline from scratch
knots = pp.knots;
coefs = pp.coefs;
% orders along all dimensions are inferred as
order1 = length(knots{1})-size(squeeze(coefs),1)
order2 = length(knots{2})-size(squeeze(coefs),2)
sizec = [order1, order2];
sp = spmak(knots,squeeze(coefs)); % somehow, supplying sizec doesn't work


figure
subplot(1,3,1)
mesh(XXe,YYe,f)
title('Approximand')
subplot(1,3,2)
mesh(XXe,YYe,fhat4)
title('Cubic')
subplot(1,3,3)
mesh(XXe,YYe,fhat4)
title('Piecewise linear')
sgtitle('spapi')

return
%% ndim_simplex
Xgrid =cell(ndim,1); % (ndim*1) cell array storing the (1*Nx) and (1*Ny) grid points
Xgrid{1} = Xq;
Xgrid{2} = Yq;
xgrid = cell(ndim,1); % (ndim*1) cell array storing the Nx*Ny nd-grid points
[xgrid{1}, xgrid{2}] = ndgrid(Xgrid{1}, Xgrid{2});

XXgrid = zeros(ndim,numel(xgrid{1})); % (ndim*Nq) double array storing the Nq interpolation nodes as (1*ndim) coordinates
for jj = 1:ndim
    XXgrid(jj,:) = xgrid{jj}(:)';
end
ns = length(XXgrid); % ns = Nq, the length of the query sample, i.e. the # of (x,y)-nodes.

x = Xgrid; % 1*ndim cell grid of query points
xx = XXgrid; % ndim * Nq  double grid of query points
ff = fq(:)'; % 1*Nq function values at query points
[alph,Q1,fracin] = ndim_simplex(x,xx,ff);
[fhat,Q2] = ndim_simplex_eval(x,xx,alph);

% post-production for plotting
fq_Ryan = reshape(fhat, size(fq));

figure
subplot(1,2,1)
mesh(XXe,YYe,f)
subplot(1,2,2)
mesh(XXq,YYq,fq_Ryan)
sgtitle('ndim\_simplex - workin y''all')

return
%% A nice little play with ndim_simplex (see Materials 33)
% Do an initial approx of the anchoring function to initialize the coeffs
ng = 10;
% grids for k^(-1)_{t-1} and f_{t|t-1}
k1grid = linspace(0.001,param.gbar,ng);
fegrid = linspace(-5,5,ng);
% values for k^{-1}_t for the grid
k = zeros(ng,ng);
for i=1:ng
    for j=1:ng
        k(i,j) = fk_smooth_pi_only(param,fegrid(j), 1./k1grid(i));
    end
end
k1 = 1./k;
% map to ndim_simplex
x = cell(2,1);
x{1} = k1grid;
x{2} = fegrid;
[xxgrid, yygrid] = meshgrid(k1grid,fegrid);
alph = ndim_simplex(x,[xxgrid(:)';yygrid(:)'],k1);
% just a check for use
xx = [k1grid(1); fegrid(1)];
out = ndim_simplex_eval(x,xx,alph);
if out==k1(1)
    disp('bingo')
end


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

% % let's check that at least this is more or less correct
% mesh(XXq,YYq,w) % yup
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Step 4. Compute Chebyshev coefficients
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
T = chebyshev(z,n);

ai=zeros(m,m);
wTT =zeros(m,m);
T2 = T.^2;
for i=0:n
    for j=0:n
        for k=1:m
            for l=1:m
                wTT(k,l) = w(k,l)*T(k,i+1)*T(l,j+1);
            end
        end
        Ti2 = sum(T2(:,i+1));
        Tj2 = sum(T2(:,j+1));
        ai(i+1,j+1) = sum(sum(wTT))/ (Ti2*Tj2);
    end
end

% might not have to compute them painfully
% ai - T\w/T; % ain't right
ai - (T'*w*T) / ((T'*T)*(T'*T)) % this is QUITE close
ai - ((T'*T)*(T'*T))\(T'*w*T) % this is also QUITE close

% Then the interpolant is
% 1. Convert the evaluation points to [-1,1] nodes
Z1 =2*(Xe-a)./(b-a)-1;
Z2 =2*(Ye-c)./(d-c)-1;
% 2. Evaluate the Chebys there
T1 = chebyshev(Z1,n);
T2 = chebyshev(Z2,n);

% 3. Interpolant
% for i=0:n
%     for j=0:n
%         p_el(i+1,j+1,:,:)=ai(i+1,j+1)*T1(:,i+1)*T2(:,j+1)';
%     end
% end
% pxy = squeeze(sum(sum(p_el,1),2));
% was again hoping not to have to compute it painfully
% pxy - T1*ai*T2' % hoppla, both of these are correct
% pxy - T2*(T1*ai)' % it's b/c T1=T2, b/c Z1=Z2. I don't know if this is a general result.
pxy = T1*ai*T2';

figure
subplot(1,2,1)
mesh(XXe,YYe,f)
subplot(1,2,2)
mesh(XXe,YYe,pxy)
sgtitle('2D Chebyshev approximation - it''s working!')



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