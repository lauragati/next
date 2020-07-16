function grid = uneven_grid(minval,maxval,ng)
% 15 July 2020

% check that ng is odd
if mod(ng,2)==0
    error('Use an odd number of breakpoints')
end

% concentrate on the positive half (the negative will mirror it)
ngpos = (ng-1)/2-1; % amount of points to place on (0,maxval)
posdist = maxval;
pospoints_evenly = 0+posdist/(ngpos+1) :posdist/(ngpos+1): maxval-posdist/(ngpos+1) ; % position of points if it was even
% Now shift the uneven points away from their even position
decreasing_shift = sort((pospoints_evenly/maxval).^(1/4), 'descend');
pospoints = abs(pospoints_evenly -decreasing_shift);
if length(pospoints)>1 && pospoints(1)>pospoints(2)
    pospoints(1) = pospoints(2)/3;
end


negpoints = sort(-pospoints);
grid = [minval, negpoints,0,pospoints, maxval];