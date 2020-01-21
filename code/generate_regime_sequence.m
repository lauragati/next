function r = generate_regime_sequence(p11,p22, T)

this_code = mfilename;
max_no_inputs = nargin(this_code);
if nargin ==0
    T = 100;
    p11 = 0.95; %0.95 (Davig and Leeper 2007 values)
    p22 = 0.93; % 0.93
elseif nargin < max_no_inputs
    T = 100;
end

p21 = 1-p11;
p12 = 1-p22;

r = nan(1,T);
r(1) = 1; % initialize in the active regime
for t=2:T
    if r(t-1)==1 % active regime
    r(t) = randsample([1,2],1,'true',[p11, p21]);
    elseif r(t-1)==2 % passive regime
    r(t) = randsample([1,2],1,'true',[p12, p22]);
    end
end