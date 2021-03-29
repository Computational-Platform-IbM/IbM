function [bac_x, bac_y] = bac_shovingloops(bac_x, bac_y, bac_r, bac_m, s_dist, overlap)
bac_n = length(bac_x);
%%% Shoving until a stable system is achieved
shov = 1;       % boolean test for more shoving (1: shoving needed; 0: shoving ready)
xr = zeros(bac_n,1);
yr = zeros(bac_n,1);
xr0 = xr; yr0 = yr;
while shov == 1 
    %%% One shoving step
    shov = 0;   % start assuming there is no need for shoving
    for i=1:bac_n-1  % browse n-1 cells
        for j=i+1:bac_n       % search overlapping cells among all other cells, "facing forward" ...
            dxx = bac_x(i,1) - bac_x(j,1) + 1e-20;
            if abs(dxx) < s_dist
                dyy = bac_y(i,1) - bac_y(j,1) + 1e-20;
                if abs(dyy) < s_dist
                    d = sqrt(dxx*dxx+dyy*dyy + 1e-20); % current distances between cell centers
                    % calculate overlap r0
                    r0 = (bac_r(i,1) + bac_r(j,1)) - d;
                    if r0 > overlap
                        shov = 1;
                        r0d = r0/d;
                        dix = (dxx/2) * r0d;
                        diy = (dyy/2) * r0d;
                        a1 = 1-bac_m(i)/(bac_m(i)+bac_m(j)) ;
                        a2 = 1-bac_m(j)/(bac_m(i)+bac_m(j));
                        xr(i) = xr(i) - dix*a1;
                        yr(i) = yr(i) - diy*a1;
                        xr(j) = xr(j) + dix*a2;
                        yr(j) = yr(j) + diy*a2;
                    end %if r0
                end %if
            end
        end %for j
    end %for i
    bac_x = bac_x - xr;
    bac_y = bac_y - yr;
    xr = xr0;
    yr = yr0;
end %while
end