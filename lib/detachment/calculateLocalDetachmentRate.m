function Fdet = calculateLocalDetachmentRate(i, j, kDet, grid, xcenter, ycenter)
    % Calculate the local detachment speed: Fdet = kDet * d^2.
    %
    % i, j: gridcell indices along x and y direction respectively
    % kDet: Detachment constant determining how fast detachment is
    % grid: struct containing all information regarding spacial
    %   discretization
    % xcenter, ycenter: x and y coordinate of the center of the granule
    %
    % -> Fdet: speed of detachment

    
    % calculate gridcell position
    x = (i - 0.5)*grid.dx;
    y = (j - 0.5)*grid.dy;
    
    % calculate distance from granule-center
    d = sqrt((x - xcenter)^2 + (y - ycenter)^2);
    
    % calculate Fdet
    Fdet = kDet * d^2;
    
    % correct to infinity if center granule is reached
    if d == 0
        Fdet = Inf;
    end
end
