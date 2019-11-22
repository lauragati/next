function fev = obj_minFEV_Ryan(gbar, xsim,ysim, gx, hx,T)
% reenact a cleaned-up sim_learn.m where I'm not letting the observables be
% determined by expectations
b = gx*hx;
b1 = b(1,:);
k = gbar^(-1);
pibar = 0; % initialize at where the simulation started (RE)
evening_fcst = nan(T,1);
morning_fcst = nan(T,1);
FEt_1 = nan(T,1);
for t = 2:T-1 % start at 2 to be compatible with the simulation                
        % Create forecasts, FE and do the updating
        morning_fcst(t) = pibar + b1*xsim(:,t); % this morning's one-step ahead forecast of tomorrow's state E(pi_{t+1} | I_{t}^m)
        FEt_1(t) = ysim(1,t)-(pibar + b1*xsim(:,t-1)); % yesterday evening's forecast error, realized today
        pibar = pibar + k.^(-1).*(ysim(1,t)-(pibar + b1*xsim(:,t-1)) );
        evening_fcst(t) = pibar + b1*xsim(:,t); % today's evening's one-step ahead forecast of tomorrow's state E(pi_{t+1} | I_{t}^e)
end
fev = nanmean(FEt_1.^2);