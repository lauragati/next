% create a grid for VFI for the midsimple anchoring model with 4 states
% 6 April 2020
% n_x = 4; the number of states
function [grid, ngrid] = create_grid_4states(N)

% Specify min and max values for each state
pibarmax = 10; pibarmin =-10;
k1max = 1; k1min= 0;
s_max = 3; s_min = -3; % gonna use same bounds for the exog states

stepsize = (pibarmax-pibarmin)/(N-1);
pgrid = pibarmin:stepsize:pibarmax;

stepsize = (k1max-k1min)/(N-1);
k1grid = k1min:stepsize:k1max;

% these are AR1s so should be discretized differently - ignore that for a
% moment
stepsize = (s_max-s_min)/(N-1);
sgrid = s_min:stepsize:s_max;


% a fake move
grid = zeros(4, N^4);
ngrid = size(grid,2);

% brute force
index=1;
for i=1:N
    pibar = pgrid(i);
    for j=1:N
        k1 = k1grid(j);
        for k=1:N
            rn = sgrid(k);
            for l=1:N
                u = sgrid(l);
                grid(:,index) = [pibar,k1,rn,u]';
                index = index+1;
            end
        end
    end
end
         

