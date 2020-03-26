function loss = objective_ACstructure(varp,param,eN,T,N,burnin,PLM,gain, acf_data, W)

this_code = mfilename;
max_no_inputs = nargin(this_code);
if nargin < max_no_inputs
    W = eye(length(acf_data));
end

% Replace the estimated parameters with the new guess
param.('d') = varp(1);
param.('c') = varp(2);

% Simulate data given the new params:
y = fun_sim_anchoring(param,T,N, burnin,eN,PLM,gain);
ny = size(y,1);

% Filter the newly simulated data

% % 1) HP filter
% g = nan(size(y));
% c = nan(size(y));
% for i=1:ny
%     [g(i,:),c(i,:)] = HPfilter(y(i,:)');
% end

% 2) Hamilton filter
% h=8;
% v = nan(ny,T-4-h+1);
% for i=1:ny
%     [v(i,:)] = Hamiltonfilter(y(i,:)');
% end

% 3) BK filter
K=12;
ystar = nan(ny,T-2*K);
for i=1:ny
    ystar(i,:) = BKfilter(y(i,:)');
end

filt=ystar;
% Compute AC structure of the newly simulated, filtered data
K=4;
acf = zeros(ny,K+1);
for i=1:ny
    for j=1:K+1
        k=j-1;
        VC = cov(filt(i,k+1:end), filt(i,1:end-k));
        acf(i,j) =VC(2,1);
    end
end
acf = vec(acf);

loss = (acf_data -acf)'*W^(-1)*(acf_data - acf);