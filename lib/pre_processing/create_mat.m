filename = './planning/Excels/Templates/AOBNOBAMXCMX_template.xlsx';

fprintf('>>>>>>>>>> LOADING EXCEL FILE\n')
[grid, bac_init, constants, settings, init_params] = loadPresetFile(filename);

fprintf('>>>>>>>>>> INITIALISING BACTERIA\n')
% initialise bacteria
n = bac_init.start_nBac; 

switch bac_init.method
    case 'granule'
        % ----- initiate circle with bacteria (random)
        radius = bac_init.granule_radius;
        tic
        [bac.x, bac.y] = blue_noise_circle(n, grid.nX/2*grid.dx, grid.nY/2*grid.dy, radius);
        toc
    case 'suspension'
        % ----- initiate suspension (semi-random)
        margin = 0.1*grid.dx*grid.nX;
        xRange = [margin, grid.dx*grid.nX - margin];
        yRange = xRange; % assume square domain
        [bac.x, bac.y] = distribute_rectangle(n, xRange, yRange);
end

bac.species = randi(length(constants.speciesNames), size(bac.x)); % random for now
bac.molarMass = 0.6 * constants.max_bac_mass_grams * ones(n, 1) / constants.bac_MW;
bac.radius = ((bac.molarMass * constants.bac_MW / constants.bac_rho) * (3 / (4 * pi))).^(1/3);
bac.active = ones(size(bac.x), 'logical');


if strcmp(bac_init.method, 'granule')
    for gg = 1:5
        bac = bacteria_shove(bac, grid, constants); % otherwise bacteria might overlap at start...
    end


    % remove cells outside of granule radius
    remove = sqrt((bac.x - grid.dx*grid.nX/2).^2 + (bac.y - grid.dy*grid.nY/2).^2) > bac_init.granule_radius;
    fprintf('%d bacteria removed outside of starting granule\n', sum(remove))
    bac.x(remove) = [];
    bac.y(remove) = [];
    bac.radius(remove) = [];
    bac.molarMass(remove) = [];
    bac.species(remove) = [];
    bac.active(remove) = [];
end

fprintf('%d starting bacteria in the system\n', length(bac.x))


% DEBUG
f = plotBacs(grid, bac, constants, 0);
% END DEBUG


function [X, Y] = rand_circle(N, x_center, y_center, r)
    Ns = round(4/pi * N + 2.5*sqrt(N) + 100);
    X = rand(Ns,1)*(2*r) - r;
    Y = rand(Ns,1)*(2*r) - r;
    I = find(sqrt(X.^2 + Y.^2)<=r);
    X = X(I(1:N)) + x_center;
    Y = Y(I(1:N)) + y_center;
end

function [X, Y] = blue_noise_circle(n, x_center, y_center, r)

    m = 20;
    X = zeros(n, 1);
    Y = zeros(n, 1);
    
    % pick initial point in center
    X(1) = x_center;
    Y(1) = y_center;
    nPoints = 1;
    
    for q = 1:n-1
        % generate set of points
        [candidate_x, candidate_y] = rand_circle(m, x_center, y_center, r);
        
        % calc distance
        dist_sq = sqrt((X(1:nPoints) - candidate_x').^2 + (Y(1:nPoints) - candidate_y').^2);% + (candidate_x'.^2 + candidate_y'.^2 - r^2);
        
        % calculate minimum distance per candidate point
        [~, i] = max(min(dist_sq, [], 1));
        
        % pick the furthest from any point
        X(nPoints+1) = candidate_x(i);
        Y(nPoints+1) = candidate_y(i);
        nPoints = nPoints + 1;        
    end

end


function [x, y] = distribute_rectangle(n, xRange, yRange)
    nSections = ceil(sqrt(n)*1.2);
    xlist = linspace(xRange(1), xRange(2), nSections)';
    ylist = linspace(yRange(1), yRange(2), nSections)';
    
    space_margin = xlist(2) - xlist(1);
    
    x = repmat(xlist, nSections, 1);
    y = kron(ylist, ones(nSections, 1));
    
    % add some random noise
    noise_x = rand(size(x)) * 0.6*space_margin - 0.3*space_margin;
    noise_y = rand(size(y)) * 0.6*space_margin - 0.3*space_margin;
    
    x = x + noise_x;
    y = y + noise_y;
    
    
    % reduce the number of bacteria
    i = randperm(length(x));
    i = i(1:n);
    
    x = x(i);
    y = y(i);
end

