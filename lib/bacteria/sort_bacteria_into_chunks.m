function bac = sort_bacteria_into_chunks(bac, grid, chunks, focus_region, nChunks_dir)
    % Reorganize the indices in the bacterial struct for easy/consecutive
    % access in parallel computing -> bacteria in same chunk next to each
    % other
    %
    % bac: struct containing all information regarding the bacteria
    % chunks: struct with start & end indices per chunk
    %
    % -> bac: see above

    % calculate which gridcell each bacterium is in
    ix = floor(bac.x / grid.dx) + 1; % +1 because of MATLAB indexing
    iy = floor(bac.y / grid.dy) + 1;

    % calculate which chunk each bacterium is in
    ixChunk = floor((ix - focus_region.x0) / chunks.dx_chunk) + 1;
    iyChunk = floor((iy - focus_region.y0) / chunks.dy_chunk) + 1;
    bac_chunk = (iyChunk - 1) + nChunks_dir * (ixChunk - 1) + 1; % why is Matlab indexing this stupid...?!
    
    % create sorting index
    [~, sortChunkIndex] = sort(bac_chunk);

    % reorganise bac struct
    bac.x = bac.x(sortChunkIndex);
    bac.y = bac.y(sortChunkIndex);
    bac.species = bac.species(sortChunkIndex);
    bac.molarMass = bac.molarMass(sortChunkIndex);
    bac.radius = bac.radius(sortChunkIndex);
    bac.active = bac.active(sortChunkIndex);
    bac.mu = bac.mu(sortChunkIndex);
end