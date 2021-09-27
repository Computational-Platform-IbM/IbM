function chunks = create_chunks(nChunks_dir, focus_region)
    % Create chunks in the focus region
    %
    % nChunks_dir: number of chunks per direction
    % focus_region: (rectangular) region of interest in the entire domain 
    %
    % chunks: struct containing indices in x & y direction per chunk

    dx = focus_region.x1 - focus_region.x0 + 1;
    dx_chunk = ceil(dx/nChunks_dir);
    
    dy = focus_region.y1 - focus_region.y0 + 1;
    dy_chunk = ceil(dy/nChunks_dir);

    % calculate which start & end indices per chunk in x direction
    indices_x = zeros(nChunks_dir, 2);
    temp_x = focus_region.x0;
    for ixChunk = 1:nChunks_dir
        indices_x(ixChunk,:) = [temp_x, temp_x + dx_chunk - 1];
        temp_x = temp_x + dx_chunk;
    end
    indices_x(end) = focus_region.x1;

    % calculate which start & end indices per chunk in x direction
    indices_y = zeros(nChunks_dir, 2);
    temp_y = focus_region.y0;
    for iyChunk = 1:nChunks_dir
        indices_y(iyChunk,:) = [temp_y, temp_y + dy_chunk - 1];
        temp_y = temp_y + dy_chunk;
    end
    indices_y(end) = focus_region.y1;
    
    % recalibrate so that indices start at 1
    %   -> assumes that focus region is already being extracted and that
    %   chunks are only additional
    indices_x = indices_x - focus_region.x0 + 1;
    indices_y = indices_y - focus_region.y0 + 1;

    chunks = struct;
    chunks.indices_x = indices_x;
    chunks.indices_y = indices_y;
    chunks.dx_chunk = dx_chunk;
    chunks.dy_chunk = dy_chunk;
end