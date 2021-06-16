
% global variables
k = 0.3; c = 35;


% create continuous space with values
coord_span = [-50, 50]; % assume square space
coord_space = linspace(coord_span(1), coord_span(2));
[X,Y] = meshgrid(coord_space);
R = sqrt(X.^2 + Y.^2);
Z = @(R) 1 ./ (1 + exp(-k*(R - c))) + eps;
f = figure(1);
f.Position = [-1773, 164, 1614, 607]; clf;
subplot(1,2,1);
pcolor(X,Y,Z(R)); 
shading interp; 
axis tight; axis equal; axis off;
colormap(viridis());
colorbar();

subplot(1,2,2);
pcolor(X,Y,Z(R)); 
shading interp; 
axis tight; axis equal; axis off;
colormap(viridis());
colorbar();
hold on;
rectangle('Position', [coord_span(1), coord_span(1), diff(coord_span), diff(coord_span)], 'FaceColor', [1,1,1,0.5], 'EdgeColor', 'none');
% create coarse grid over the space & assign values to grid
nX = 10; dx = diff(coord_span) / nX;
nY = 10; dy = diff(coord_span) / nY;
grid_values_level0 = zeros(nX, nY);
grid_centers_level0 = zeros(nX, nY, 2);

for xi = 1:nX
    for yi = 1:nY
        x0 = coord_span(1) + dx * (xi - 1);
        y0 = coord_span(1) + dy * (yi - 1);
        w = dx;
        h = dy;
                
        % assign value
        x_center = x0 + w/2;
        y_center = y0 + h/2;
        grid_centers_level0(xi, yi, :) = [x_center, y_center];
        r = sqrt(x_center^2 + y_center^2);
        grid_values_level0(xi, yi) = Z(r);
        
        % draw rectangle for grid cell in figure
        rectangle('Position', [x0, y0, w, h], 'EdgeColor', 'k', 'LineWidth',3);
    end
end

% automagically divide too coarse cells into finer cells (once)
flag_for_division_0 = zeros(nX, nY);
for xi = 1:nX
    for yi = 1:nY
                
        % up
        yii_up = min(yi + 1, nY); xii_up = xi;
                
        % down
        yii_down = max(yi - 1, 1); xii_down = xi;
        
        % right
        yii_right = yi; xii_right = min(xi + 1, nX);
        
        % left
        yii_left = yi; xii_left = max(xi - 1, 1);
        
        flag_for_division_0(xi, yi) = needsDivision(grid_values_level0, xi, yi, xii_up, yii_up) || ...
            needsDivision(grid_values_level0, xi, yi, xii_down, yii_down) || ...
            needsDivision(grid_values_level0, xi, yi, xii_left, yii_left) || ...
            needsDivision(grid_values_level0, xi, yi, xii_right, yii_right);
    end
end

% calculate values in level1 and draw gridcells
divM = [-1, 1; 1, 1; 1, -1; -1, -1];
grid_values_level1 = zeros(nX, nY, 4);
grid_centers_level1 = zeros(nX, nY, 4, 2);
for xi = 1:nX
    for yi = 1:nY
        if flag_for_division_0(xi, yi)
            centers_new = reshape(grid_centers_level0(xi, yi, :), [1,2]) + [dx, dy]./4 .* divM;
            grid_centers_level1(xi, yi, :, :) = centers_new;
            r = sqrt(sum(centers_new.^2, 2));
            grid_values_level1(xi, yi, :) = Z(r);

            for q = 1:4
                % draw rectangle for grid cell in figure
                rectangle('Position', [centers_new(q,1)-w/4, centers_new(q,2)-h/4, w/2, h/2], 'EdgeColor', 'k', 'LineWidth',2);
            end
        end
    end
end

% make it recursive, so that cells can then divide again
% automagically divide too coarse cells into finer cells (once)
flag_for_division_1 = zeros(nX, nY, 4);
for xi = 1:nX
    for yi = 1:nY
        if flag_for_division_0(xi, yi)
            for q = 1:4
                self = grid_values_level1(xi, yi, q);
                % for now: assume that boundary cells are not divided, thus
                % not checking with min/max
                if q == 1 % nw
                    % up
                    yii = yi + 1; xii = xi; qq = 4;
                    val_up = grid_values_level1(xii, yii, qq);
                    % down
                    yii = yi; xii = xi; qq = 4;
                    val_down = grid_values_level1(xii, yii, qq);
                    % right
                    yii = yi; xii = xi; qq = 2;
                    val_right = grid_values_level1(xii, yii, qq);
                    % left
                    yii = yi; xii = xi - 1; qq = 2;
                    val_left = grid_values_level1(xii, yii, qq);
                elseif q==2 % ne
                    % up
                    yii = yi + 1; xii = xi; qq = 3;
                    val_up = grid_values_level1(xii, yii, qq);
                    % down
                    yii = yi; xii = xi; qq = 3;
                    val_down = grid_values_level1(xii, yii, qq);
                    % right
                    yii = yi; xii = xi + 1; qq = 1;
                    val_right = grid_values_level1(xii, yii, qq);
                    % left
                    yii = yi; xii = xi; qq = 1;
                    val_left = grid_values_level1(xii, yii, qq);
                elseif q==3 % se
                    % up
                    yii = yi; xii = xi; qq = 2;
                    val_up = grid_values_level1(xii, yii, qq);
                    % down
                    yii = yi - 1; xii = xi; qq = 2;
                    val_down = grid_values_level1(xii, yii, qq);
                    % right
                    yii = yi; xii = xi + 1; qq = 4;
                    val_right = grid_values_level1(xii, yii, qq);
                    % left
                    yii = yi; xii = xi; qq = 4;
                    val_left = grid_values_level1(xii, yii, qq);
                elseif q==4 % se
                    % up
                    yii = yi; xii = xi; qq = 1;
                    val_up = grid_values_level1(xii, yii, qq);
                    % down
                    yii = yi - 1; xii = xi; qq = 1;
                    val_down = grid_values_level1(xii, yii, qq);
                    % right
                    yii = yi; xii = xi; qq = 3;
                    val_right = grid_values_level1(xii, yii, qq);
                    % left
                    yii = yi; xii = xi - 1; qq = 3;
                    val_left = grid_values_level1(xii, yii, qq);
                end
                
                flag_for_division_1(xi, yi, q) = needsDivisionLvl(self, val_up) || ...
                    needsDivisionLvl(self, val_down) || ...
                    needsDivisionLvl(self, val_right) || ...
                    needsDivisionLvl(self, val_left);
            end
        end
    end
end

grid_values_level2 = zeros(nX, nY, 4, 4);
grid_centers_level2 = zeros(nX, nY, 4, 4, 2);
for xi = 1:nX
    for yi = 1:nY
        for q = 1:4
            if flag_for_division_1(xi, yi, q)
                centers_new = reshape(grid_centers_level1(xi, yi, q, :), [1,2]) + [dx, dy]./8 .* divM;
                grid_centers_level2(xi, yi, q, :, :) = centers_new;
                r = sqrt(sum(centers_new.^2, 2));
                grid_values_level2(xi, yi, q, :) = Z(r);

                for qq = 1:4
                    % draw rectangle for grid cell in figure
                    rectangle('Position', [centers_new(qq,1)-w/8, centers_new(qq,2)-h/8, w/4, h/4], 'EdgeColor', 'k', 'LineWidth',1);
                end
            end
        end
    end
end


% plot 1 grid solution (original code)
figure();
pcolor(X,Y,Z(R)); 
shading interp; 
axis tight; axis equal; axis off;
colormap(viridis());
colorbar();
hold on;
rectangle('Position', [coord_span(1), coord_span(1), diff(coord_span), diff(coord_span)], 'FaceColor', [1,1,1,0.5], 'EdgeColor', 'none');
% create coarse grid over the space & assign values to grid
nX = 40; dx = diff(coord_span) / nX;
nY = 40; dy = diff(coord_span) / nY;

for xi = 1:nX
    for yi = 1:nY
        x0 = coord_span(1) + dx * (xi - 1);
        y0 = coord_span(1) + dy * (yi - 1);
        w = dx;
        h = dy;
                
        % draw rectangle for grid cell in figure
        rectangle('Position', [x0, y0, w, h], 'EdgeColor', 'k', 'LineWidth',1);
    end
end



function bool = needsDivision(vals, xi, yi, xii, yii)
%%%% function to determine whether a division is required based on two
%%%% neighbouring nodes
    cutoff = 0.7;

    neighb = vals(xii, yii);
    self = vals(xi, yi);

    bool = abs(neighb - self) / (min(neighb, self) + 0.1) > cutoff;
end

function bool = needsDivisionLvl(self, neighb)
%%%% function to determine whether a division is required based on two
%%%% neighbouring nodes
    cutoff = 0.7;
    if neighb == 0
        bool = false;
    else
        bool = abs(neighb - self) / (min(neighb, self) + 0.1) > cutoff;
    end
end