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

chebyshev = @(Tn,Tn_1,n,x) 2.*x.*Tn -Tn_1;

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
            varphi(k,j) = legendre(varphi(k-1,j), varphi(k-2,j),k,x(j));
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
% This completes the quadrature.

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
            varphi(k,i) = legendre(varphi(k-1,i), varphi(k-2,i),k,x_i);
        end
    end
end
px = bet_ols'*varphi;

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% 2.) Chebyshev interpolation

% Note: number of nodes (m) = number of Chebyshev coefficients (n+1)
% Step 1: Compute ngrid+1 interpolation nodes
k=1:N+1;
z_k = cos((2.*k - 1)/(2*N+2) * pi);
% Step 2. Adjust nodes to [a,b] interval (which is idiotic here but I'm doing it for completeness)
x_k = (z_k + 1).*(b-a)/2 + a;
% Step 3: evaluate f at the nodes
y_k = true_fun(x_k);
% Step 4: compute Chebyshev coefficients a_i
a_i = zeros(N+1,1);
Ti = zeros(N+1,N+1);
% First loop to construct Chebyshev coefficients
for i=1:N+1
    % first evaluate the chebyshev
    if i==1 % order 0
        Ti(i,:) = 1;
    elseif i==2 % order 1
%         Ti(i,:) = x_k;
        Ti(i,:) = z_k;
    else % order 2 to N-1
        Ti(i,:) = chebyshev(Ti(i-1,:), Ti(i-2,:),i,z_k);
    end
    sum_num = 0;
    sum_den = 0;
    for k=1:N+1
        sum_num = sum_num + y_k(k)*Ti(i,k);
        sum_den = sum_den + Ti(i,k)^2;
    end
    a_i(i) = sum_num/sum_den;
end
% Pass through all the points in the interval to create the big Chebyshev
Ti = zeros(N+1,T);
for j=1:T
    x_j = interval(j);
    % Evaluate Chebyshev of order k at that node x(j) recursively
    for i=1:N+1
        if i==1 % order zero
            Ti(i,j) = 1;
        elseif i==2 % order 1
            Ti(i,j) = 2.*(x_j-a)/(b-a)-1;
        else % order 2 to N
            Ti(i,j) = chebyshev(Ti(i-1,j), Ti(i-2,j),1,2.*(x_j-a)/(b-a)-1);
        end
    end
end

fhat =0;
for i=1:N+1
    fhat = fhat + a_i(i) .* Ti(i,:);
end

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% 3.) Cubic spline
% We have ngrid interpolation points, ngrid-1 intervals, so 4(ngrid-1)
% coefficients. Let's first determine the nodes.
node_indexes = 1:T/(ngrid-1):T;
node_indexes(end+1) = T; % add the ending node (note: they are not quite evenly spaced)
xi = interval(node_indexes);
% Evaluate f at the nodes
yi = true_fun(xi);
% Now initial guess for the coefficients
coeffs0 = zeros(ngrid-1,4); % [a,b,c,d]

% set up objective spline with equation system
resids = obj_spline_secant_hermite(coeffs0,xi,yi);

% %Optimization Parameters
options = optimoptions('fsolve', 'TolFun', 1e-9, 'display', 'iter', 'MaxFunEvals', 10000);


tic
objh = @(coeffs) obj_spline_secant_hermite(coeffs,xi,yi);
[coeffs_opt,FVAL] = fsolve(objh,coeffs0, options);
toc

% fish out coefficients
a = coeffs_opt(:,1);
b = coeffs_opt(:,2);
c = coeffs_opt(:,3);
d = coeffs_opt(:,4);


% construct spline
sx = zeros(T,1);
segment_index=1;
for i=1:T
    x=interval(i); 
    if i>=node_indexes(segment_index+1) && segment_index < ngrid-1
        segment_index = segment_index+1;
    end
    sx(i) = a(segment_index)+ b(segment_index)*x + c(segment_index)*x^2 + d(segment_index)*x^3;
end
    

% Plot configs
[fs, lw] = plot_configs;

figure
set(gcf,'color','w'); % sets white background color
set(gcf, 'Position', get(0, 'Screensize')); % sets the figure fullscreen
h1 = plot(interval, f,'linewidth',lw); hold on
h2 = plot(interval, px, '--','linewidth',lw);
h3 = plot(interval, fhat, ':','linewidth',lw);
h4 = plot(interval, sx,'linewidth',lw);
ax = gca; % current axes
ax.FontSize = fs;
grid on
grid minor
title('True f against approximations', 'FontSize',fs)
legend([h1,h2, h3, h4], 'Approximand', 'Legendre least squares', 'Chebyshev interpolation', 'Spline', 'FontSize',fs, 'location', 'northwest')
legend('boxoff')


if isequal(func2str(true_fun), '@(x)(x+1).^(1/4)')
    figname='example2';
else
    figname='example3';
end
if print_figs ==1
    disp(figname)
    cd(figpath)
    export_fig(figname)
    cd(current_dir)
    close
end


