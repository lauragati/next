% anchoring_functional_forms.m
% just experiment with functional forms for g
rho=0.9;
gam=0.001;
delt=0.01;
fe = randn(T,1);
k(1) = 0.145;
for t=2:T
k(t) = rho*k(t-1) +gam*(fe(t)^2);
end

plot(k)