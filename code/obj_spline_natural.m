function eqs = obj_spline_natural(coeffs,x,y)

% Unpack the coefficients
a = coeffs(:,1);
b = coeffs(:,2);
c = coeffs(:,3);
d = coeffs(:,4);
ngrid = size(coeffs,1)+1;
n=ngrid-1; % number of intervals

% Now set up the equation system
eqs = zeros(4*n,1);
index=0;

% 1.) Interpolating conditions + continuity at interior nodes 
% the first conditions (actually this is Judd, 6.9.2)
for i=1:n
    index=index+1;
eqs(index) = -y(i) + a(i)+ b(i)*x(i) + c(i)*x(i)^2 + d(i)*x(i)^3;
end
% the second conditions (actually this is Judd, 6.9.1)
for i=2:n+1
    index=index+1;
eqs(index) = -y(i) + a(i-1)+ b(i-1)*x(i) + c(i-1)*x(i)^2 + d(i-1)*x(i)^3;
end

% check the number is right
if index ~= 2*n
    warning('wrong number of equations')
end

% 2.) Twice differentiability at interior nodes
for i=2:n
    index=index+1;
    % first derivatives
    eqs(index) = -b(i) - 2*c(i)*x(i) -3*d(i)*x(i)^2    +b(i-1) + 2*c(i-1)*x(i) +3*d(i-1)*x(i)^2; 
    index=index+1;
    % second derivatives
    eqs(index) = - 2*c(i) - 6*d(i)*x(i)    + 2*c(i-1) + 6*d(i)*x(i-1); 
end

% natural spline last 2 conditions
index = index+1;
eqs(index) = -b(1)-2*c(1)*x(1)-3*d(1)*x(1)^2;
% to simplify, let's call 
index=index+1;
eqs(index) = -b(n)-2*c(n)*x(n+1)-3*d(n)*x(n+1)^2;