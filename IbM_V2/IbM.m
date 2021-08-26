function IbM(options)
%%%%% main function to run the IbM model with given (optional) options 
% output_directory: folder in which to store the results and temporary
%   variables. Creates it if not present already.
% preset: if given a preset file, create the model
%   from that, otherwise checks in the output_directory for an R.mat file to
%   load from.
% path_shoving: path to the java (.jar) file used for shoving

    %% argument validation
    arguments
        options.output_directory {checkFolder} = 'Results';
        options.preset {mustBeFile, mustBeFormat(options.preset, {'xls', 'xlsx', 'mat'})}
        options.path_shoving {mustBeFile, mustBeFormat(options.path_shoving, {'jar'})} = [pwd '\lib\shovingQuadTree.jar']; 
    end
    

    %% java import for shoving
    javaaddpath(options.path_shoving);
    addpath(genpath('lib')); % make every subfolder with functions accessible to the code


    
    %% preset loading
    R = loadPreset(options);
    
    
    %% enable/disable debug disp/warning
    warning('on', 'DEBUG:noActionRequired');
    warning('on', 'DEBUG:actionRequired');
    
    
    %% ========== Time advancements ==========
    fprintf('> SIMULATION RUNNING >>>>>\n');
    
    R = integTime(R, options.output_directory); %#ok<NASGU>
    
    fprintf('> SIMULATION FINISHED >>>>>\n');
    

end

function R = loadPreset(options)
    % if preset is given, then load from preset (either process from excel
    % or load R.mat)
    if isfield(options, 'preset')
        if isFormat(options.preset, {'mat'})
            R = load(options.preset);
        else
            R = loadModelXlsx(options.preset);
        end

    % if no preset is given, then check in the output directory for a R.mat
    % file and load from that else throw error
    else
        if exist([options.output_directory, '\R.mat'], 'file')
            R = load([options.output_directory, '\R.mat']);
        else
            eid = 'Preset:notFound';
            msg = 'No preset file was found. Either supply one in the function call, or add one in the output directory.';
            throwAsCaller(MException(eid, msg))
        end
    end
end

function checkFolder(dir)
    if ~isfolder(dir)
        mkdir(dir);
    end
end

function mustBeFormat(filename, formats)
    if ~isFormat(filename, formats)
        eid = 'File:formatIncorrect';
        msg = ['File <', strrep(filename, '\', '\\'), '> must be of the format: .', char(join(formats, ', .'))];
        throwAsCaller(MException(eid, msg))
    end
end

function bool = isFormat(filename, formats)
    filename_split = split(filename, '.');
    bool = any(strcmp(formats, filename_split(end)));
end

