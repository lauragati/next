function B = projection_approx(seq,n_input_jumps,param,gx,hx,eta,PLM,gain,T,ndrop,e, alph_opt,x, k1grid,fegrid, g_fe)
% Compute steps 3 and 4 of the projection procedure.
% 16 May 2020
kapp = param.kapp;
sig  = param.sig;
psi_pi = param.psi_pi;
psi_x  = param.psi_x;
rho_k = param.rho_k;
gam_k = param.gam_k;
lamx = param.lamx;
alph= param.alph;
bet = param.bet;

[s1, s2, s3, s4] = smat(param,hx);

% Step 2 - compute synthetic time series v
% Simulate given input sequences
% input k
% [xsim, ysim, k, ~, FA, FB, FEt_1] = sim_learnLH_clean_given_seq2(param,gx,hx,eta,PLM, gain, T,ndrop,e,seq,n_input_jumps);
% input FE
[xsim, ysim, k, phi_seq, FA, FB, FEt_1,g_pi,g_pibar] = sim_learnLH_clean_given_seq3_approx(param,gx,hx,eta,PLM, gain, T,ndrop,e,seq,n_input_jumps,...
    alph_opt,x, k1grid,fegrid, g_fe);
pi = ysim(1,2:end-1);
x  = ysim(2,2:end-1);
xl = ysim(2,1:end-2);
i  = ysim(3,2:end-1);
s  = xsim(:,2:end-1);
st_1 = xsim(:,1:end-2);
fa = FA(:,2:end-1);
fb = FB(:,2:end-1);
fe = FEt_1(1,2:end-1); % select the FE(pi);
k1 = 1./k(2:end-1);
pibar = squeeze(phi_seq(1,1,1:end-2))';
b = gx*hx;
g_pi = g_pi(2:end-1);
g_pibar = g_pibar(2:end-1);


% Create basis as 1st, 2nd and 3rd powers of the states
X = [k1;pibar;s([1,end],:)];
sx = [X', X'.^2, X'.^3];

% % Step 3. Compute analogue synthetic expectation series
% tau=T-2;
% E = zeros(tau,1);
% for t=1:tau
%     for i=1:(tau-t)
%         PI=1;
%         for j=0:i-1
%             PI = PI*(1 - k1(t+1+j) - fe(t+1+j)*g_pibar(t+j));
%         end
%         SIG = x(t+i)*PI;
%     end
%     E(t) = SIG;
% end

% Step 3. Alternative: following Ryan's comment of 22 May 2020, use means
% of variables to form expectations (in order to avoid introducing trends due to truncation)
meanx  = mean(x);
meank1 = mean(k1);
meanfe = mean(fe);
meangb = mean(g_pibar);
tau=T-2;
E = zeros(tau,1);
for t=1:tau
    for i=1:(tau-t)
        PI=1;
        for j=0:i-1
            PI = PI*(1 - meank1 - meanfe*meangb);
        end
        SIG = meanx*PI;
    end
    E(t) = SIG;
end

% Step 4.
B = (sx'*sx)\(sx'*E);










