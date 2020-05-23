function [tv, pibar, k1] = mat31_TV(OM,param,pp,i,pibart_1,k1t_1, s, st_1,sgrid,PI)
Om6  = OM{6};
Om7  = OM{7};
Om9  = OM{9};
Om10 = OM{10};
Om11 = OM{11};
bet =param.bet;
ns = length(sgrid);

% period loss
Lt =  Om6*i^2 + Om7*pibart_1^2 + Om9*i*s + Om10*pibart_1*s + Om11*s;

% Implied endogenous states at t
k1 = mat31_k1(OM,param,k1t_1,i,pibart_1,s,st_1);
pibar = mat31_pibar(param,OM,pibart_1,k1,s,st_1,i);

% approximation of value function % here form weighted average
% expectations over tomorrow's r and u
v = zeros(ns,ns);
for i=1:ns
    rp=sgrid(i);
    for j=1:ns
        up = sgrid(j);
        v(i,j) = fnval(pp,{k1,pibar,rp,up});
    end
end

tv = Lt + bet*(v(:)'*PI(:));

