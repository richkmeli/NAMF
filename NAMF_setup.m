% NAMF_setup.m
% Configuration script for NAMF - Numerical Analysis MATLAB Functions
% Automatically adds all necessary paths to the MATLAB path

function NAMF_setup()
    % Get the current directory path
    current_dir = fileparts(mfilename('fullpath'));
    
    % Add all subdirectories to the path
    addpath(genpath(current_dir));
    
    % Confirmation message
    fprintf('NAMF: All paths have been added to the MATLAB path\n');
    fprintf('Available functions:\n');
    fprintf('  - Quadrature formulas: %d functions\n', count_m_files('Quadrature_Formulas'));
    fprintf('  - Root finding: %d functions\n', count_m_files('Root_Finding'));
    fprintf('  - Linear systems: %d functions\n', count_m_files('Linear_Systems'));
    
    % Save the path for future sessions
    savepath;
    fprintf('Path saved for future sessions.\n');
end

function count = count_m_files(folder)
    % Count .m files in a folder and subfolders
    files = dir(fullfile(folder, '**/*.m'));
    count = length(files);
end