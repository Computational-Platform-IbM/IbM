function plotDetachTime(grid, T)
    % plot the detachment times per gridcell.
    %
    % grid: struct containing all information regarding space
    %   discretization
    % T: matrix (nX*nY) with the time of detachment per grid cell
    

    figure(13);
    
    try % try updating first, otherwise draw figure entirely
        ax = gca;
        ax.Children(1).CData = T';
        caxis([0, 12]);
        drawnow();
    catch e
        switch e.identifier
            case {'MATLAB:noSuchMethodOrField', 'MATLAB:noPublicFieldForClass'} % has not been made before
                imagesc([grid.dx/2, grid.nX*grid.dx - grid.dx/2], [grid.dy/2, grid.nY*grid.dy - grid.dy/2], T'); hold on;
                colormap(viridis());
                colorbar();
                caxis([0, 12])
                axis equal;
                xlim([0, grid.nX*grid.dx]);
                ylim([0, grid.nY*grid.dy]);
                set(gca, 'YDir', 'normal'); 
                title('Detachment times [h]');
            otherwise
                rethrow(e)
        end
    end
end

