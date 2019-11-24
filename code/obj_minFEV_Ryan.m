function fev = obj_minFEV_Ryan(gbar, xsim,ysim, gx, hx,T)
% reenact a cleaned-up sim_learn.m where I'm not letting the observables be
% determined by expectations
b = gx*hx;
b1 = b(1,:);
k = gbar^(-1);
pibar = 0; % initialize at where the simulation started (RE)
FEt_1 = nan(T,1);
for t = 2:T % start at 2 to be compatible with the simulation                
        % Create forecasts, FE and do the updating
        FEt_1(t) = ysim(1,t)-(pibar + b1*xsim(:,t-1)); % yesterday evening's forecast error, realized today
        pibar = pibar + k^(-1)*  (ysim(1,t)-(pibar + b1*xsim(:,t-1)) );
end
fev = nanmean(FEt_1.^2);