% materials27
% do Judd, Numerical approximation examples 2-3 on 233 Mac.

clearvars
close all
clc

% Add all the relevant paths and grab the codename
this_code = mfilename;
[current_dir, basepath, BC_researchpath,toolpath,export_figpath,figpath,tablepath,datapath] = add_paths;
todays_date = strrep(datestr(today), '-','_');

% Variable stuff ---
print_figs        = 0;
stop_before_plots = 0;
skip_old_plots    = 0;
output_table = print_figs;

skip = 1;


%% Let's compute some polynomials (Judd, Numerical, p. 211 Mac)
tau = 400;
interval = linspace(-1,1,tau);
N = 10; % degree of polynomial (0 to N-1)

legendre = @(Pn,Pn_1,n,x) (2*n+1)/(n+1)*x*Pn - n/(n+1)*Pn_1;

chebyshev = @(Tn,Tn_1,n,x) 2*x*Tn -Tn_1;

laguerre = @(Ln,Ln_1,n,x) 1/(n+1)*(2*n+1-x)*Ln -n/(n+1)*Ln_1;

hermite = @(Hn,Hn_1,n,x) 2*x*Hn -2*n*Hn_1;

P = zeros(N,tau);
T = zeros(N,tau);
L = zeros(N,tau);
H = zeros(N,tau);

for i=1:tau
    x=interval(i);
    for n=0:N-1
        
        if n==0 % order 0
            P(n+1,i) =1;
            T(n+1,i) =1;
            L(n+1,i) =1;
            H(n+1,i) =1;
        elseif n==1 % order 1
            P(n+1,i) = x;
            T(n+1,i) = x;
            L(n+1,i) = 1-x;
            H(n+1,i) = 2*x;
        else % order 2 to N-1
            P(n+1,i) = legendre(P(n,i), P(n-1,i),n,x);
            T(n+1,i) = chebyshev(T(n,i), T(n-1,i),n,x);
            L(n+1,i) = laguerre(L(n,i), L(n-1,i),n,x);
            H(n+1,i) = hermite(H(n,i), H(n-1,i),n,x);
        end
    end
end


figure
subplot(2,2,1)
plot(interval,P(2,:)); hold on
plot(interval,P(3,:))
plot(interval,P(4,:))
plot(interval,P(5,:))
title('Legendre')
subplot(2,2,2)
plot(interval,T(2,:)); hold on
plot(interval,T(3,:))
plot(interval,T(4,:))
plot(interval,T(5,:))
title('Chebyshev')
subplot(2,2,3)
plot(interval,L(2,:)); hold on
plot(interval,L(3,:))
plot(interval,L(4,:))
plot(interval,L(5,:))
title('Laguerre')
subplot(2,2,4)
plot(interval,H(2,:)); hold on
plot(interval,H(3,:))
plot(interval,H(4,:))
plot(interval,H(5,:))
title('Hermite')

% yay they work!



%% Example 2: approximand = (x+1)^(1/4) and  Example 3: approximand = min(max(-1,4*(x-0.2),1))
close all
N=4;
ngrid =5; % number of points used in approximation
T =100; % "size of dataset"
interval = linspace(-1,1,T);
a  = interval(1);
b  = interval(end);

% Approximand in "dataset"
% true_fun = @(x) (x+1).^(1/4);
% To do Example 3 uncomment the following
true_fun = @(x) min(max(-1*ones(size(x)),4.*(x-0.2)),1*ones(size(x)));
f = true_fun(interval);

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% 1.) Legendre least squares
% thank god the Legendre weights are just 1
om_l = @(x) 1;

x = zeros(ngrid,1);
omq = zeros(ngrid,1);
varphi=zeros(N+1,ngrid); % order N plus 1 because of order 0
truth = zeros(ngrid,1);
element_num = zeros(ngrid,1);
element_den = zeros(ngrid,1);
bet_ols =zeros(N+1,1);
% Gauss-Chebyshev quadrature to evaluate the inner products
% First loop: evaluate all the quadrature nodes, weights, Legendre
% polynomials of all orders up to N, and approximand
% For ngrid nodes
for j=1:ngrid
    % Calculate node
    x(j) = a + (b-a)*(cos((2*j-1)/(2*ngrid)*pi)+1)/2;
    % Calculate quadrature weights
    omq(j) = (b-a)/2 * pi/ngrid * sqrt(1-cos((2*j-1)/(2*ngrid)*pi)^2);
    
    % Evaluate Legendre of order k at that node x(j) recursively
    for k=1:N+1
        if k==1 % order zero
            varphi(k,j) = 1;
        elseif k==2 % order 1
            varphi(k,j) = x(j);
        else % order 2 to N
            varphi(k,j) = legendre(varphi(k-1,j), varphi(k-2,j),1,x(j));
        end
    end
    
    % approximand at xj
    truth(j) = true_fun(x(j));
end

% Second loop: evaluate inner products for each order k to create beta_ols
for k=1:N+1
    for j=1:ngrid
        % Calculate each element of the inner products for order k
        element_num(j) = omq(j) * truth(j)*varphi(k,j)*om_l(x(j));
        element_den(j) = omq(j) * varphi(k,j)*varphi(k,j)*om_l(x(j));
    end
    inner_prod_k_num = sum(element_num);
    inner_prod_k_den = sum(element_den);
    bet_ols(k) = inner_prod_k_num/inner_prod_k_den;
end

% Third loop: evaluate the Legendre least squares as the fitted value over
% the "dataset". First create Legendre polynomials over the dataset
varphi = zeros(k,T);
for i=1:T
    x_i = interval(i);
    % Evaluate Legendre of order k at that node x(j) recursively
    for k=1:N+1
        if k==1 % order zero
            varphi(k,i) = 1; 
        elseif k==2 % order 1
            varphi(k,i) = x_i;
        else % order 2 to N
            varphi(k,i) = legendre(varphi(k-1,i), varphi(k-2,i),1,x_i);
        end
    end
end
px = bet_ols'*varphi;

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% 2.) Chebyshev interpolation

figure
h1 = plot(interval, f); hold on
h2 = plot(interval, px, '--');
h3 = plot(interval, zeros(size(interval)), ':');
h4 = plot(interval, zeros(size(interval)));
title('True f against approximations')
legend([h1,h2, h3, h4], 'Approximand', 'Legendre least squares', 'Chebyshev interpolation', 'Spline')




