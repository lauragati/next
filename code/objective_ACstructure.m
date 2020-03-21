function loss = objective_ACstructure(varp,param,eN,T,N,burnin,PLM,gain, acf_data, W)

this_code = mfilename;
max_no_inputs = nargin(this_code);
if nargin < max_no_inputs
    W = eye(length(acf_data));
end

% Replace the estimated parameters with the new guess
param.('d') = varp;

% Simulate data given the new params:
y = fun_sim_anchoring(param,T,N, burnin,eN,PLM,gain);
ny = size(y,1);

% Filter the newly simulated data
K=12;
ystar = nan(ny,T-2*K-1);
for i=1:ny
    ystar(i,:) = BKfilter(y(i,:)');
end

% Compute AC structure of the newly simulated, filtered data
K=4;
acf = zeros(ny,K+1);
for i=1:ny
    for j=1:K+1
        k=j-1;
        VC = cov(ystar(i,k+1:end), ystar(i,1:end-k));
        acf(i,j) =VC(2,1);
    end
end
acf = vec(acf);

loss = (acf_data -acf)'*W^(-1)*(acf_data - acf);