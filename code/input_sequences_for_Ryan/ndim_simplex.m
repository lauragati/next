% NDIM_SIMPLEX - Compute the linear finite element approximation of a function
% using grid points x, samples points xx, and values f.
%
% usage
%
% [alph,Q] = ndim_simplex(x,xx,f)
%
% where
%
% x = 1-by-N cellarray, each elemnt constaining a the nx(i) grid points on
% dimension i, for i = 1...N
% xx = an N-by-ns set of sample points for the function
% f  = an 1-by-ns set of sample values for the function
%
% alph,Q are such that f = Q*alph if ns = prof(nx). Otherwise, Q*alph
% minimizes sum of squared residuals, sum((f-Q*alph).^2).
%
% Ryan Chahrour, Boston College, 2014, ryan.chahrour@bc.edu.
%
% generalized from Weiser and Zarantonello, "A note on piecewise linear and multilinear table interpolation
% in many dimensions" (1988)

function [alph,Q,fracin] = ndim_simplex(x,xx,f)

%Dimensions
N = length(x);          %Number of dimensions
nx = length(xx);        %Total number of evaluated points
nzq = (N+1)*nx;         %Non-zero entries in Q matrix

%Initializations
idx = zeros(size(xx));  %Bins for each value on each grid
gg = zeros(size(xx));   %Normalized coordinates within bin
nn = zeros(1,N);        %Number of grid points for each dimenion

Qvals = zeros(1,nzq);           %Values for Q matrix
Qidx  = zeros(1,nzq);           %Column indexes for Q matrix
xidx = repmat((1:nx)',[1,N+1]); %Row indexes for Q matrix

for jj = 1:N
    %Select grid points
    x_ex = x{jj};
    
    %Binning to handle values beyond grid...
    xtmp      = x_ex;
    xtmp(1)   = -inf;
    xtmp(end) = inf;
    
    %number of point on each grid
    nn(jj) = length(x_ex);
    
    %Bin all values
    [~,idx(jj,:)] = histc(xx(jj,:),xtmp);
    
    %Local coordinates for each point
    gg(jj,:) = (xx(jj,:) - x_ex(idx(jj,:)))./(x_ex(idx(jj,:)+1)-x_ex(idx(jj,:)));
end
   

%Sorting each point
[gg,xord] = sort(gg,1);


%Carefully vectorized!
idx = idx';
gg  = gg';
xord = xord';


%Computing the weights on each element
nncum = repmat([1,cumprod(nn(1:end-1))],[nx,1]);

keepx = ((max(gg)<=1) + (min(gg)>=0)) == 2;
fracin = sum(keepx)/length(keepx);
%This section is designed to find weights for each node
%using the definition of S in Weiser-Zarantonello 1988.
%The weights are x(1), (x(j+1) - x(j)), 1-x(N+1)


%First entry is is x(1)
S = ones(nx,N);
idx0 = (sum(nncum.*(idx-1+S),2))+1;
Qidx(1:nx) = idx0;
Qvals(1:nx) = gg(:,1);

%For middle grids, x(j+1)-x(j)
for jj = 1:N-1
    S(nx*(xord(:,jj)-1)+(1:nx)') = 0;  %Increment down S matrix
    idx1 = (sum(nncum.*(idx-1+S),2))+1;
    Qidx((nx*jj+1):nx*(jj+1)) = idx1;
    Qvals((nx*jj+1):nx*(jj+1)) = gg(:,jj+1)-gg(:,jj);
end

%For final, 1-x(N+1)
idxn = (sum(nncum.*(idx-1),2))+1;
Qidx(N*nx+1:end) = idxn;
Qvals(N*nx+1:end) = 1-gg(:,end);

%Create sparse matrix
Q = sparse(xidx(:)',Qidx,Qvals,nx,prod(nn),nzq);

%Evaluate best-fit weights
alph = Q\f(:);