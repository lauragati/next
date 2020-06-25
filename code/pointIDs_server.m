ngrid = 6;
pointID = zeros(6,ngrid^ngrid);
tic
datestr(now)
iter=0;
for i=1:ngrid
    for j=1:ngrid
        for k=1:ngrid
            for l=1:ngrid
                for m=1:ngrid
                    for n=1:ngrid
                        iter=iter+1;
                        pointID(:,iter) = [i;j;k;l;m;n];
                    end
                end
            end
        end
    end
end
toc
save('6by6loss_pointIDs.mat',pointID)