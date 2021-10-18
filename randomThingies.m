% testing speed of extracting data from conc matrix
% squeeze vs reshape
% row first or column first

%{

Important take-home message:
matlab is column-major, meaning that the last index changes slowest (first
index values lay next to each other in memory).

so to get the concentration of the entire domain (diffusion) it is fastest
to have x,y,c -> conc = conc(:,:,ic)

and for the extraction of all compound concentration in one grid cell it is
fastest to have c,x,y -> conc = conc(ic, :, :)

MAIN POINT: when looping in column-major -> loop from last index to first
index i.e. in (i, j, k) loop as:
    for k
        for j
            for i
                do stuff...
            end
        end
    end

%}



clc;

x = 1024;
y = x;
c = 8;
conc_xyc = rand(x, y, c);

tic
for i = 1:x
    for j = 1:y
        a = squeeze(conc_xyc(i, j, :));
    end
end
timeSqueeze = toc;

tic
for i = 1:x
    for j = 1:y
        a = reshape(conc_xyc(i, j, :), [], 1, 1);
    end
end
timeReshape_xy = toc;

tic
for j = 1:x
    for i = 1:y
        a = reshape(conc_xyc(i, j, :), [], 1, 1);
    end
end
timeReshape_yx = toc;

fprintf('\nConc: xyc\n')
fprintf('Extracting concentration at (x, y) \n\tSqueeze: %.4f s\n\tReshape (xy): %.4f s\n\tReshape (yx): %.4f s\n', timeSqueeze, timeReshape_xy, timeReshape_yx)

% optimize data storage in memory for data extraction
conc_cyx = permute(conc_xyc, [3, 2, 1]);

tic
for i = 1:x
    for j = 1:y
        a = squeeze(conc_cyx(:, i, j));
    end
end
timeSqueeze = toc;

tic
for i = 1:x
    for j = 1:y
        a = reshape(conc_cyx(:, i, j), [], 1, 1);
    end
end
timeReshape_xy = toc;

tic
for j = 1:y
    for i = 1:x
        a = conc_cyx(:, i, j);
    end
end
timeNative_yx = toc;

tic
for i = 1:x
    for j = 1:y
        a = conc_cyx(:, i, j);
    end
end
timeNative = toc;

fprintf('\nConc: cyx\n')
fprintf('Extracting concentration at (x, y) \n\tSqueeze: %.4f s\n\tReshape: %.4f s\n\tNative (xy):  %.4f s\n\tNative (yx):  %.4f s\n', timeSqueeze, timeReshape_xy, timeNative, timeNative_yx)



tic
for ci = 1:c
    comp_conc = reshape(conc_cyx(ci, :, :), x, y);
end
timeDiff_cyx = toc;

tic
for ci = 1:c
    comp_conc = conc_xyc(:, :, ci);
end
timeDiff_xyc = toc;

fprintf('\nDiffusion slicing\n')
fprintf('Extracting concentration over entire domain \n\tXYC: %.4f s\n\tCYX: %.4f s\n',timeDiff_xyc, timeDiff_cyx)

