% materials27
% do Judd, Numerical approximation examples on 233 Mac.

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

interval = linspace(-1,1,10);
N_x = numel(interval);
N = 5; % degree of polynomial
ngrid = 5; % number of points used in approximation
%% Let's compute some polynomials (Judd, Numerical, p. 211 Mac)
legendre = @(Pn,Pn_1,n,x) (2*n+1)/(n+1)*x*Pn - n/(n+1)*Pn_1;

chebyshev = @(Tn,Tn_1,n,x) 2*x*Tn -Tn_1;

laguerre = @(Ln,Ln_1,n,x) 1/(n+1)*(2*n+1-x)*Ln -n/(n+1)*Ln_1;

hermite = @(Hn,Hn_1,n,x) 2*x*Hn -2*n*Hn_1;

P = zeros(N,N_x);
T = zeros(N,N_x);
L = zeros(N,N_x);
H = zeros(N,N_x);

for i=1:N_x
    x=interval(i);
    for n=0:N-1
        
        if n==0
            P(n+1,i) =1;
            T(n+1,i) =1;
            L(n+1,i) =1;
            H(n+1,i) =1;
        elseif n==1
            P(n+1,i) = x;
            T(n+1,i) = x;
            L(n+1,i) = 1-x;
            H(n+1,i) = 2*x;
        else
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

%% Unfortunately, least squares approx requires the computation of inner products, which are integrals. So first implement some numerical integration
% Judd, Numerical, p 258 Mac.
% Do  it for ff=x^(1/4)
ff = @(x) x^(1/4);
% ff = @(x) exp(x);
simpson = @(a,b) (b-a)/6 * (ff(a) + 4*ff((a+b)/2) + ff(b));
a=0;
b=1;
% the simple simpson rule is very crude
simpson(a,b);

% so now do a composite rule which breaks the interval into n smaller
% subintervals of length h
n=14; % n+1 points
h = (b-a)/n;
x0 = a;
f0 = ff(x0);
clear x
clear f
clear summand
for j=1:n
    x(j) = a+j*h;
    % Evaluate approximand at all nodes
    f(j) = ff(x(j));
    % check with what it should be multiplied
    if j==n
        m=1;
    elseif mod(j,2)==0
        m=2;
    elseif mod(j,2)==1
        m=4;
    end
    summand(j) = m*f(j);
end
composite_simpson =h/3 *(f0+ sum(summand))

% you know, I'm not exactly reproducing the table in Judd, but honestly I think he
% makes mistakes with what is n, whether it's the number of points or the
% number of points minus 1.

return
clear f
clear x

%% Example 1: approximand = e^(2x+2)

true1 = @(x) exp(2*x + 2);
f1 = true1(interval);

% 1.) Legendre least squares
% thank god the Legendre weights are just 1
om = @(x) 1;

% For each polynomial order k
for k=1:N
    % each inner product is an integral
    % use Gauss-Chebyshev quadrature (see tryouts.m)
    a  = interval(1);
    b  = interval(end);
    n  = N;
    % Then for n nodes
    for j=1:n
        % Calculate node
        x(j) = a + (b-1)*(cos((2*j-1)/(2*n)*pi)+1)/2;
        % Evaluate Legendre of order k at that node x(j) recursively
        varphi=zeros(k,1);
        for l=1:k
            if l==1
                varphi(l) =1;
            elseif l==2
                varphi(l) = x(j);
            else
                varphi(l) = legendre(varphi(l-1), varphi(l-2),1,x(j));
            end
        end
        % Calculate quadrature weights
        omega(j) = (b-a)/2 * pi/n * sqrt(1-cos((2*j-1)/(2*n)*pi)^2);
        
        % Calculate each element of the inner products
        element_num(j) = omega(j)*true1(x(j))*varphi(end)*om(x(j));
        element_den(j) = omega(j)*varphi(end)*varphi(end)*om(x(j));
    end
    inner_product_num = sum(element_num);
    inner_product_den = sum(element_den);
    OLS = inner_product_num/inner_product_den * varphi;
end
px = sum(OLS);





return
figure
plot(interval, f1)
title('True f against approximations')


%% Example 2: approximand = (x+1)^(1/4)

true2 = @(x) (x+1).^(1/4);
f2 = true2(interval);

figure
plot(interval,f2)
title('True f against approximations')



%% Example 3: approximand = min(max(-1,4*(x-0.2),1))

true3 = @(x) min(max(-1*ones(size(x)),4.*(x-0.2)),1*ones(size(x)));
f3 = true3(interval);

figure
plot(interval,f3)
title('True f against approximations')