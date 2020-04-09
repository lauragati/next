% evaluate all residuals in the model using an RE-discretion targeting
% criterion
% 8 April 2020
function resids = res_RE_TC(param,hx,ysim,k,s, fa,fb)
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

T = length(k);
pi = ysim(1,:);
x  = ysim(2,:);
i  = ysim(3,:);


resTC = zeros(T-2,1);
resA9 = zeros(T-2,1);
resA10 = zeros(T-2,1);

for t=2:T-1
   % 1.) Evaluate RE-discretion TC residuals
    resTC(t) = -pi(t) - lamx/kapp *x(t);
    % 2.) Evaluate IS and PC residuals
    resA9(t) = -x(t) -sig*i(t) + stuff1*fb(:,t) + stuff2*s(:,t);
    resA10(t) =-pi(t) +kapp*x(t) + stuff3*fa(:,t) + stuff4*s(:,t);
end

% resTC = -pi -lamx/kapp*x;
% resA9 =  -x -sig*ysim(3,:) + stuff1*fb + stuff2*s;
% resA10 =-pi +kapp*x + stuff3*fa + stuff4;



resids = [resA9'; resA10'; resTC'];