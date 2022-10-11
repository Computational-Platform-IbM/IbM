clear all;
filename = './planning/Excels/main.xlsx';
javaaddpath([pwd '/lib/bacteria/shovingQuadTreekDist.jar']);

fprintf('>>>>>>>>>> LOADING EXCEL FILE\n')
[grid, bac_init, constants, settings, init_params] = loadPresetFile(filename);

molarMass = 0.6 * constants.max_bac_mass_grams / constants.bac_MW;
radius = ((molarMass * constants.bac_MW / constants.bac_rho) * (3 / (4 * pi))).^(1/3);

fprintf('>>>>>>>>>> INITIALISING BACTERIA\n')
% initialise bacteria

switch settings.model_type
    case {'granule', 'mature granule'}
        n = bac_init.start_nBac; 

        % ----- initiate circle with bacteria (random)
        radius_granule = bac_init.granule_radius;
        tic
        [bac.x, bac.y] = blue_noise_circle(n, grid.nX/2*grid.dx, grid.nY/2*grid.dy, radius_granule);
        toc
    case 'suspension'
%         n = bac_init.start_nColonies * bac_init.start_nBacPerColony; 

        % ----- initiate suspension (semi-random)
        margin = 0.2*grid.dx*grid.nX; % 20% of simulation domain as margin for growth (empirical)
        xRange = [margin, grid.dx*grid.nX - margin];
        yRange = xRange; % assume square domain
        
        r_colony = sqrt((bac_init.start_nBacPerColony * radius * constants.kDist)^2) / 5; % empirical... 
        [bac.x, bac.y] = distribute_microcolonies(bac_init.start_nColonies, bac_init.start_nBacPerColony, r_colony, xRange, yRange);
end

bac.molarMass = ones(numel(bac.x), 1) * molarMass;
bac.radius = ones(numel(bac.x), 1) * radius; 
bac.active = ones(size(bac.x), 'logical');

for gg = 1:5
    bac = bacteria_shove(bac, grid, constants); % otherwise bacteria might overlap at start...
end

if ismember(settings.model_type, {'granule', 'mature granule'})
    % remove cells outside of granule radius
    remove = sqrt((bac.x - grid.dx*grid.nX/2).^2 + (bac.y - grid.dy*grid.nY/2).^2) > bac_init.granule_radius;
    fprintf('%d bacteria removed outside of starting granule\n', sum(remove))
    bac.x(remove) = [];
    bac.y(remove) = [];
    bac.radius(remove) = [];
    bac.molarMass(remove) = [];
%     bac.species(remove) = []; % assigned after shoving (due to mature
%     granule possibility...)
    bac.active(remove) = [];
end

if strcmp(settings.model_type, 'mature granule')
    bac.species = AMXinside(bac, grid, constants);
else
    bac.species = randi(length(constants.speciesNames), size(bac.x)); % random for now
end

fprintf('%d starting bacteria in the system:\n', length(bac.x))
for s = 1:numel(constants.speciesNames)
    fprintf('\t%d %s\n', sum(bac.species == s), constants.speciesNames{s})
end

fprintf('>>>>>>>>>> DONE!\n')


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


function [x, y] = distribute_microcolonies(nColonies, nBacPerCol, r_colony, xRange, yRange)
    nSections = ceil(sqrt(nColonies)*1.1);
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
    i = i(1:nColonies);
    
    x = x(i);
    y = y(i);
    
    % create microcolony at each x, y
    for ix = 1:numel(x)
        [x_microcol, y_microcol] = blue_noise_circle(nBacPerCol, x(ix), y(ix), r_colony);
        x = [x; x_microcol(2:end)];
        y = [y; y_microcol(2:end)];
    end
end


function species = AMXinside(bac, grid, constants)
    % calculate distance from center of granule
    dist_center = sqrt((bac.x - grid.nX/2*grid.dx).^2 + (bac.y - grid.nY/2*grid.dy).^2);
    [~, I] = sort(dist_center);
    
    % randomly select how many AMX in the system
    nAMX = sum(randi(length(constants.speciesNames), size(bac.x)) == 3);
    species = zeros(size(bac.x));
    
    % make individuals closest to center of the domain AMX
    iAMX = find(strcmp(constants.speciesNames, 'AMX')); % AMX -> species = ?
    species(I(1:nAMX)) = iAMX;
    
    % assign all other species
    species_other = randi(length(constants.speciesNames) - 1, size(species(species == 0)));
    species_other(species_other >= iAMX) = species_other(species_other >= iAMX) + 1;
    species(I(nAMX+1:end)) = species_other;
end
