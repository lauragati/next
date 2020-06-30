function [param, set, param_names, param_values_str, param_titles] = parameters_next
param.bet  = 0.99;%0.99
param.sig  = 1; %1 IES=1 (log utility, consistent with balanced growth) Peter: micro estimates= 1/4 or 1/5, asset pricing lit finds 1/100.
param.alph = 0.9; %0.9,0.5 (prob that firm stuck with price), set to correspond to an expected duration of 2 quarters. Collard: estimates of alpha between 0.66-0.75 (adjust prices every 3-4 months)
param.eta = 1/4; %1/4, inverse of Frisch. won't matter here. The Frisch: 1 or 2 (micro) or 4 (macro) e.g. Basu Lec 9, slide 16 Mac.
param.om = 1.25;% 1.25 Subsumes all things that go into elasticity of marginal cost to output (including Frisch) Woodford. Interest and Prices, p. 172 (Table 3.1).
param.thet = 10; %10. price elasticity of demand, Woodford, taken from Chari Kehoe McGrattan 2000.
param.zeta = (param.om + 1/param.sig)/(1+param.om*param.thet); % parameter of strategic complementarity. If < 1, strat comp in price setting. If >1, strt subs.
param.kapp = param.zeta * (1-param.alph*param.bet)/param.alph; % Woodford. Interest and Prices, p. 187.
param.psi_x  = 0; %0.01,0, 1
param.psi_pi = 1.5; %1.5
param.w = 1+param.sig*param.psi_x +param.kapp*param.sig*param.psi_pi;
param.gbar    = 0.145; % 0.145 param_correct CEMP. 0.02 is the value a dgain algorithm gets after 50 periods.
param.thetbar = 16;%16, 1 or 4 or 0.029; % param_correct CEMP
param.rho_r = 0; %0
param.rho_i = 0.6; % 0.6 CEMP: 0.877, not a perfect mapping, but the MC process, standing in for demand. Too much: set 0.7 to get monpol shock to increase i on impact
param.rho_u = 0; %0
param.rho = 0; % persistence of lag of interest rate
param.sig_r = 1;%1 to facilitate IRFs. 0.1; %?
param.sig_i = 1; %0 to turn it off for optimal monpol. 1 to facilitate IRFs.  0.359 = sig_e from CEMP, standing in for the demand shock
param.sig_u = 1; %1 to facilitate IRFs.  0.277 = sig_mu from CEMP, the cost-push shock
param.kap =  0.8; % 0.8 allows the CUSUM test-statistic to be revised at a different rate than the estimate of the mean.  0 < kap < 1.
param.thettilde =2.5;%1.6 or 2.5 the new thetbar for CUSUM-test. I just set both to match CEMP's criterion approx.
% param.gam = 0.128; %0.128 indexation in NKPC. Posterior mean CEMP.
% param.p11 = 0.95; %0.95 Follow Davig and Leeper's transition probabilities (Prob active|active = 0.95, Prob(passive|passive)=0.93)
% param.p22 = 0.93;%0.93
% param.p21 = 1-param.p11;
% param.p12 = 1-param.p22;
% param.psi1 = 1.8; % 2.19 Taylor-coefficient on inflation in acctive regime (Davig and Leeper 2007 values)
% param.psi2 = 0.89; % 0.89 Taylor-coefficient on inflation in acctive regime (Davig and Leeper 2007 values)
param.lamx = 0;%0.01, % 0 Rotemberg Woodford 1997 estimate 0.05. Woodford 2011 suggests optimal value = kapp/theta (0.01683)
param.lami = 0;
% param.d = 10; % 10 slope of anchoring function.
% param.c = 0; % 0 intercept of anchoring function.
param.rho_k = 0.6151; %0.5, 0.9 or 0.5 Anchoring function parameters for Alternative 7.
param.gam_k = 0.0011; %0.001 or 0.01. Estimates LOM gain (co): 0.6151    0.0011, estimates constant only, pi only (1, 0.0047)
% param.psi_k     = 0.01; % reaction function coeffs for a first-pass reaction function r1.
% param.psi_pibar = 0.01;
% param.psi_xbar  = 0.01;


% % "weighted il" extension
% param.psi_x  = param.psi_x*param.rho;
% param.psi_pi = param.psi_pi*param.rho;

% Set is just a copy of param as of 16 Nov 2019
param_values = struct2array(param);
param_values_str = cell(size(param_values));
param_names = cell(size(param_values));
fn = fieldnames(param);
for k=1:numel(fn)
    if( isnumeric(param.(fn{k})) )
        % do stuff
        set.(fn{k}) = param.(fn{k}); % assign set the same value as param
        %         fn{k} % this prints the name of the field
        param_names{k} = fn{k};
        param_values_str{k} = replace(num2str(param_values(k)), '.','_');
    end
end

param_titles = {'\beta', '\sigma', '\alpha','\eta', '\omega', '\theta','\zeta','\kappa','\psi_x', '\psi_{\pi}',...
    'w', '\bar{g}', '\bar{\theta}','\rho_r', '\rho_i', '\rho_u','\rho','\sigma_r', '\sigma_i',...
    '\sigma_u', '\tilde{kappa}', '\tilde{\theta}','\gamma', 'p_{11}', 'p_{22}','p_{21}','p_{12}', '\psi_{1}',...
    '\psi_2','\alpha^{CB}'};
% % create indexes for the positions of the parameters
% make_index(fn)

