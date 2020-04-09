% evaluate all residuals in the midsimple anchoring model
% 3 April 2020
% updated 4 April 2020
% checked thru 9 April 2020
function resids = res_anchTC(param,hx,ysim,k,phi_seq,s,g_pi, fa,fb, H)
lamx = param.lamx;
kapp = param.kapp;
bet  = param.bet;
alph = param.alph;
sig  = param.sig;
nx = length(hx);

stuff1 = [sig,1-bet,-sig*bet]; % fa(3) I hope is rational
% stuff1 = [sig,1-bet,0];
stuff2 = sig*[1,0,0]*(eye(nx)-bet*hx)^(-1);
stuff3 = [(1-alph)*bet,kapp*alph*bet,0];
stuff4 = [0,0,1]*(eye(nx)-alph*bet*hx)^(-1);

T = length(k)-H;
pi = ysim(1,:);
x  = ysim(2,:);
pibar = squeeze(phi_seq(1,1,:));
b = squeeze(phi_seq(:,2:end,1));

resTC = zeros(T,1);
resA9 = zeros(T,1);
resA10 = zeros(T,1);

for t=2:T
    
    % 1.) evaluate target criterion residual
    sumprod = x(t+1); % when i=1
    for i=2:H-1
        prod = 1;
        for j=1:i-1
            prod = prod*(1-k(t+1+j)^(-1)-(pi(t+1+j)-pibar(t+j)-b(1,:)*s(:,t+j))*g_pi(t+j)); % this assumed perfect foresight - not if I input the expected shocks in simulation
            %             prod = prod*(1-k(t+1+j)^(-1)-(pi(t+1+j)-pibar(t+j) -b(1,:)*hx^(j-1)*s(:,t))*g_pi(t+j)); % RE CB expectations for states
        end
        sumprod = sumprod + x(t+i)*prod;
    end
    
    resTC(t) = - pi(t) - lamx/kapp*x(t) ...
        +lamx/kapp*(1-alph)*bet/(1-alph*bet)*(k(t)^(-1) + (pi(t) - pibar(t-1) - b(1,:)*s(:,t-1))*g_pi(t))...
        *sumprod;
    
    % 2.) Evaluate IS and PC residuals
    resA9(t) = -x(t) -sig*ysim(3,t) + stuff1*fb(:,t) + stuff2*s(:,t);
    resA10(t) =-pi(t) +kapp*x(t) + stuff3*fa(:,t) + stuff4*s(:,t);
end


resids = [resA9'; resA10'; resTC'];