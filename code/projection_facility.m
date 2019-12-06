function keep = projection_facility(phi,q)
% q is the threshold number eigenvals need to be smaller than. Usually this
% will be 1, reflecting that the macroeconomy is stationary.

max_no_inputs = nargin('projection_facility');
if nargin < max_no_inputs %no shock specified
    q=1; % default case
end

keep =  max(max(abs(eig(phi)))) < q;
% If keep=1, the new phi is used, otherwise, it is discarded and the old
% phi is used.