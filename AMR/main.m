clear all;

% figures:
% 1 = concentration field (continuous)
% 2 = grid visualisation
% 3 = concentration field (discretised)
% 4 = grid refinement per iteration
% 5 = progress on concentration in slice through y=0


%% define simulation domain
coord_span = 15; % assume square space
grid.xmin = -coord_span; grid.xmax = coord_span;
grid.ymin = -coord_span; grid.ymax = coord_span;
grid.nX = 32; grid.nY = 32; % number of CELLS per direction


%% define AMR parameters
grid.maxLevel = 3;


%% plot concentration gradient
visualise_continuous_concentration(grid, 0);
clim = get_limit_from_figure1(); % get this from figure 1


%% initialize grid
[grid, nCells, nNodes] = initialise_grid(grid);

% visualise_grid(grid, nCells);
% visualise_discretised_concentration(grid, nCells, clim, 0);

[grid, nCells] = refine_grid(grid, nCells, nNodes, 0); % initial refinement of grid, with exact solution


visualise_grid(grid, nCells);
visualise_discretised_concentration(grid, nCells, clim, 0);

% %% integrate over time with ADI or Crank-Nicolson (if no refinement is done)
% t = 0;
% dt = 0.5;
% old_conc = grid.cells.conc(1:nCells);
% 
% for i = 1:ceil(4/dt)
% %     grid.cells.conc(1:nCells) = ADI(grid, dt, nCells);
%     grid.cells.conc(1:nCells) = Crank_Nicolson(grid, dt, nCells);
%     
%     t = t + dt;
% %     visualise_slice_discrete(grid, t);
% end
% 
% visualise_discretised_concentration(grid, nCells, clim, t);

%% solve diffusion timestep
t = 0;
dt = 0.5;

for i = 1:ceil(4/dt)
    grid.cells.conc(1:nCells) = AMR_Diffusion(grid, nCells, dt);
    t = t + dt;
    visualise_slice_discrete(grid, t);
%     [grid, nCells] = refine_grid(grid, nCells, nNodes, t); % initial refinement of grid, with exact solution

end

visualise_discretised_concentration(grid, nCells, clim, t);


%% solve diffusion with multigrid method

% conc = multigrid_diffusion(grid, nCells, dt);


%% funtions
function visualise_grid(grid, nCells)
    f = figure(2); 
    clf;
    f.Position = [-469 556 425 400]; % fig1: [-1773 164 700 600]
    
    % find dx and dy for base level
    dy = grid.dy;
    dx = grid.dx;
    
    % draw base grid natively (no rectangles, but just lines)
    x=grid.xmin:dx:grid.xmax;
    line([x;x], [grid.ymin;grid.ymax].*ones(size(x)), 'Color', 'k', 'LineWidth', 1); hold on;
    y=grid.ymin:dy:grid.ymax;
    line([grid.xmin;grid.xmax].*ones(size(y)), [y;y], 'Color', 'k', 'LineWidth', 1);

    
    % for all higher level cells, draw them over top
    for ci = (grid.nX*grid.nY+1):nCells
        if ~grid.cells.children(ci, 1)
           % find bottom left node
           ni = grid.cells.nod(ci, 1);
           x = grid.nodes.x(ni);
           y = grid.nodes.y(ni);
           
           % draw rectangle for cell
           rectangle('Position', [x, y, dx/(2^grid.cells.lvl(ci)), dy/(2^grid.cells.lvl(ci))], 'EdgeColor', 'k', 'LineWidth', 1);
        end
    end
    
    axis equal;
    xlim([grid.xmin, grid.xmax]);
    ylim([grid.ymin, grid.ymax]);
end

function visualise_discretised_concentration(grid, nCells, clim, t)
    f = figure(3); 
    clf;
    f.Position = [-1179 86 700 840]; % fig1: [-1773 164 700 600]
    sgtitle('Discretised concentration');
    sf1 = subplot(2,1,1);
    title('2D concentration profile')
    % find dx and dy for base level
    dy = grid.dy;
    dx = grid.dx;
    
    % create color scale
    cmap = colormap(viridis());
        
    % for all cells, draw concentration
    for ci = 1:nCells
        if ~grid.cells.children(ci,1)
           % find bottom left node
           ni = grid.cells.nod(ci, 1);
           x = grid.nodes.x(ni);
           y = grid.nodes.y(ni);
           frac = (grid.cells.conc(ci) - clim(1))/diff(clim);
           c = cmap(min(max(round(length(cmap)*frac), 1), 256),:);           
           
           % draw rectangle for cell
           rectangle('Position', [x, y, dx/(2^grid.cells.lvl(ci)), dy/(2^grid.cells.lvl(ci))], 'FaceColor', c, 'EdgeColor', 'none', 'LineWidth', 0.5);
        end
    end
    
    axis equal; axis off; 
    colorbar();
    caxis(clim);
    xlim([grid.xmin, grid.xmax]);
    ylim([grid.ymin, grid.ymax]);
    
    % start at the half point of the first column of cells (assume an even
    % number of cells -> otherwise other algorithm required)
    c0 = grid.nY / 2;
    if floor(c0) ~= c0
        error('An even number of cells in each direction is required');
    end
    
    x_array = [];
    conc_array = [];
    lvl_array = [];
    
    % get first node
    ni = grid.cells.nod(c0, 2);
    x0 = grid.nodes.x(ni);
    qq = 0;
    while grid.nodes.cel(ni, 3) ~= 0
        % find the cells on either side of the y-axis
        cells = grid.nodes.cel(ni, [3,4]);
        if grid.cells.lvl(cells(1)) ~= grid.cells.lvl(cells(2))
            error('Cells on either side of the y-axis are not of the same level');
        end
        ni = grid.cells.nod(cells(1), 4);
        x1 = grid.nodes.x(ni);
        x_array(end+1) = (x1 + x0)/2;
        conc_array(end+1) = (grid.cells.conc(cells(1)) + grid.cells.conc(cells(2)))/2;
        lvl_array(end+1) = grid.cells.lvl(cells(1));
        
        x0 = x1; 

        qq = qq + 1;
    end
    
    
%     xx = linspace(grid.xmin, grid.xmax);
%     plot(xx, xx, 'LineWidth', 2);
    
    % plot concentration through center of granule
    sf2 = subplot(2,1,2);
    coord_space = linspace(grid.xmin, grid.xmax);
    plot(coord_space, GaussianDiffusionExact(coord_space, 0, t), 'k', 'LineWidth', 1); hold on;
    scatter(x_array, conc_array, [], lvl_array, 'o', 'LineWidth', 2);
    cmap = [12, 174, 85; 51, 59, 151; 236, 33, 91]./255;
    colormap(gca, cmap);
    legend('Exact solution', 'Discretised concentration');
    title(sprintf('Concentration through y=0 @ t=%.2f', t));

    
    reposition_subplots(sf1, sf2);
end

function visualise_slice_discrete(grid, t)
    f = figure(5); 
    clf;
    f.Position = [-1179 86 700 350];
    
    % start at the half point of the first column of cells (assume an even
    % number of cells -> otherwise other algorithm required)
    c0 = grid.nY / 2;
    if floor(c0) ~= c0
        error('An even number of cells in each direction is required');
    end
    
    x_array = [];
    conc_array = [];
    
    % get first node
    ni = grid.cells.nod(c0, 2);
    x0 = grid.nodes.x(ni);
    qq = 0;
    while grid.nodes.cel(ni, 3) ~= 0
        % find the cells on either side of the y-axis
        cells = grid.nodes.cel(ni, [3,4]);
        if grid.cells.lvl(cells(1)) ~= grid.cells.lvl(cells(2))
            error('Cells on either side of the y-axis are not of the same level');
        end
        ni = grid.cells.nod(cells(1), 4);
        x1 = grid.nodes.x(ni);
        x_array(end+1) = (x1 + x0)/2;
        conc_array(end+1) = (grid.cells.conc(cells(1)) + grid.cells.conc(cells(2)))/2;
        
        x0 = x1; 

        qq = qq + 1;
    end
    
    % plot concentration through center of granule
    plot(x_array, conc_array, 'o'); hold on;
    coord_space = linspace(grid.xmin, grid.xmax);
    plot(coord_space, GaussianDiffusionExact(coord_space, 0, t));
    title(sprintf('Discretised concentration along y=0 @ t=%.2f', t))
    xlabel('Position');
    ylabel('Concentration');
end


function visualise_continuous_concentration(grid, t)
    % create continuous space with values
    coord_space = linspace(grid.xmin, grid.xmax);
    [X,Y] = meshgrid(coord_space);
    
    f = figure(1);
    f.Position = [-1889 86 700 840];%[-1773, 164, 1614, 607]; 
    clf;
    sgtitle('Exact concentration');

    sf1 = subplot(2,1,1);
    pcolor(X, Y, GaussianDiffusionExact(X, Y, t)); 
    shading interp; 
    axis tight; axis equal; axis off;
    colormap(viridis());
    h = colorbar();
    h.Limits(1) = 0;
    title('2D concentration profile')


    sf2 = subplot(2,1,2);
    plot(coord_space, GaussianDiffusionExact(coord_space, 0, 0), 'LineWidth', 2); hold on;
    title('Concentration through y=0');

    reposition_subplots(sf1, sf2);
end

function conc = Z(X, Y)
    % global variables
    k = 0.3; c = 35;

    R = sqrt(X.^2 + Y.^2);
    conc = 1 ./ (1 + exp(-k*(R - c))) + eps;
end

function lim = get_limit_from_figure1()
    f1 = figure(1);
    colorBar1 = f1.Children(2);
    lim = colorBar1.Limits;
end

function reposition_subplots(sf1, sf2)
    % reposition subplots
    ht = 0.8;
    sp = (1-ht)/4; % one extra spacer for title at the top
    h1 = 3 * ht / 4 ;
    h2 = ht - h1;
    w = 0.7;
    dw = 0.12;
    pause(0.001)
    sf2.Position = [dw 1.2*sp w h2];
    sf1.Position = [dw 2.2*sp+h2 w h1];
    pause(0.001)
    sf2.Position = [dw 1.2*sp w h2];

end

