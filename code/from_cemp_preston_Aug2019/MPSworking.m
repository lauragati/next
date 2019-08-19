function [pibar_smooth, k_smooth, xi_smooth] = MPSworking(Y,param,N,M,xn_seq, xl_seq,P_seq,qtilde_seq)
% Marginalized particle smoother of CEMP for a full dataset. See their
% appendix and my notes.
% This is the cleaned-up version of marginalized_smoother_test_1by1.m
% 8 August 2019.
% Logic:
% A.) Run the forward filter, storing at each t the swarms (pibar_t, k_t, w_t xi_t) and P_tt.
% (Call them pibar_seq, k_seq, w_seq, xi_seq, P_seq).
% B.) Run the backward filter (smoother) for t=T-1:-1:1
% 1.) Calculate w_smooth
% 2a.) Draw (pibar_smooth, k_smooth) at each t from the swarm (pibar_seq,
% k_seq), using w_smooth as weights.
% 2b.) Draw xi_smooth from a normal given (pibar_smooth, k_smooth, xi_seq)


% Notation:
% N = of CEMP, number of particles. (indexed by i, or at the most, k)
% M = M of CEMP, number of histories. (indexed by j)

T = size(Y,2);
% l = size(Y,1);

% Preallocate the stored stuff from the filter
k_seq = squeeze(xn_seq(:,1,:)); % (T x N)
pibar_seq = squeeze(xn_seq(:,2,:)); % (T x N)
xi_seq = xl_seq; % (T x k_l x N)
w_seq = qtilde_seq; % (T X N)

[~, ~, ~, ~, set, variable] = parameters_cemp;
% Here I choose which parameters are estimated (which are variable, not
% constant).
% The idea is to set everything to the constant values 'set', and then
% replace the variable ones with the input parameter values.

variables = find(variable);
if isempty(variables)==1
    paramset = param;
else
    paramset = struct2array(set);
    paramset(variables) = param;
end

thetbar = paramset(2);
gbar    = paramset(3);
gam     = paramset(4);
Gam     = paramset(5);
sige    = paramset(7);
sigmu = paramset(8);

Sxi = [1,1; 1,0; 1,1];
SIG = diag([sige sigmu]);
Q = vertcat(zeros(1,4), horzcat(zeros(3,1), Sxi*SIG*Sxi'));
w_smooth = zeros(M,N);
pibar_smooth = zeros(T,M);
k_smooth = zeros(T,M);
xi_smooth = zeros(T,3,M);
% Initialize them at the end of forward filter:
initial = randsample(1:N,M,'true',w_seq(T,:));
pibar_smooth(T,:) = pibar_seq(T,initial);
k_smooth(T,:) = k_seq(T,initial);
xi_smooth(T,:,:) = squeeze(xi_seq(T,:,initial)) + squeeze(P_seq(T,:,:))^(1/2)*randn(3,M);

countwnan =0;
for t=T-1:-1:1
    for j=1:M
        % 1.) Calculate w_smooth
        sumwproba=0;
        proba=zeros(1,N);
        for i=1:N
            I1 = k_smooth(t+1,j) == gbar^(-1);
            I2 = abs((1-gam)*(Gam-1)*pibar_seq(t,i)) > thetbar*(sige+sigmu);
            I3 = abs(k_seq(t,i) - (k_smooth(t+1,j)-1)) < 1e-10 ;
            I = I1*I2 + (1-I1)*(1-I2)*I3;
            if I==1
                [~, fpibar, Apibar, fxi, Axi] = cemp_functions_matrices(param, k_seq(t,i), pibar_seq(t,i), 3);
                f = [fpibar; fxi];
                A = [Apibar; Axi];
                Omtilde = Q + A*squeeze(P_seq(t,:,:))*A';
                data = [pibar_smooth(t+1,j),xi_smooth(t+1,:,j)]';
                MU = f + A*xi_seq(t,:,i)';
                proba_ji = mvnpdf(data,MU,Omtilde);
            else
                proba_ji = 0;
            end
            proba(i)=proba_ji;
            sumwproba = sumwproba + w_seq(t,i)*proba(i);
        end
        for i=1:N
            w_smooth(j,i) = w_seq(t,i)*proba(i)/sumwproba;
        end
        % 2a.) Draw (pibar_smooth, k_smooth) at each t from the swarm (pibar_seq,
        % k_seq), using w_smooth as weights.
        index = randsample(1:N,1,'true',w_smooth(j,:));
        pibar_smooth(t,j) = pibar_seq(t,index);
        k_smooth(t,j) = k_seq(t,index);
        
        % 2b.) Draw xi_smooth from a normal given (pibar_smooth, k_smooth, xi_seq)
        [~, fpibar, Apibar, fxi, Axi] = cemp_functions_matrices(param, k_smooth(t,j), pibar_smooth(t,j), 3);
        fs = [fpibar; fxi];
        As = [Apibar; Axi];

        Deltaj = squeeze(P_seq(t,:,:))*As'* (Q + As*squeeze(P_seq(t,:,:))*As')^(-1);
        Lambdaj = squeeze(P_seq(t,:,:)) - Deltaj*As*squeeze(P_seq(t,:,:));
        normalmean = xi_seq(t,:,index)' + Deltaj*([pibar_smooth(t+1,j), xi_smooth(t+1,:,j)]' - fs - As*xi_seq(t,:,index)');
        
        try
            xi_smooth(t,:,j) = mvnrnd(normalmean,Lambdaj,1);
        catch e 
            testlam4 = Lambdaj + diag(abs(diag(Lambdaj))+0.000001 - diag(Lambdaj));
            testlam4 = (testlam4 + testlam4')/2;
            xi_smooth(t,:,j) = mvnrnd(normalmean,testlam4,1);
        end
    end
end


