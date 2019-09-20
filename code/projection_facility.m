function keep = projection_facility(phi,q)
% q is the threshold number eigenvals need to be smaller than. Usually this
% will be 1, reflecting that the macroeconomy is stationary.

keep =  max(max(abs(eig(phi)))) < q;
% If keep=1, the new phi is used, otherwise, it is discarded and the old
% phi is used.