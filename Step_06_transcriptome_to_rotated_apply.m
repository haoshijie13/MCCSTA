% Add necessary toolboxes to MATLAB path
addpath '/spm12' % SPM12 for neuroimaging processing

% Define data directories
fixfolder = '/mnt/Data18Tb/.../T402-528'; % Registered 3D RNA data
geneexpfolder = '/mnt/Data18Tb/.../segment_coronal/'; % Gene expression data
movfolder = '/mnt/Data18Tb/.../T402-528'; % Original 2D RNA slices

% Configure output paths
check_dir = './coronal_check/'; % QC directory
savepath = '/segment_coronal_to_3D'; % Final output

% File list initialization
fixfiles = dir(fullfile(fixfolder, '*.nii.gz'));
movfiles = dir(fullfile(movfolder, '*.nii.gz'));
geneexpfiles = dir(fullfile(geneexpfolder, '*.nii.gz'));

% Main processing loop
for i = 1:numel(movfiles)
    % File identification
    curfixID = fixfiles(i).name;
    curID = curfixID(1:4); % Extract slice ID (e.g. T402)
    curmovID = movfiles(i).name;
    curgeneID = geneexpfiles(i).name;

    % Load imaging data
    fiximg = squeeze(niftiread(fullfile(fixfolder, curfixID))); % 3D reference
    movimg = squeeze(niftiread(fullfile(movfolder, curmovID))); % 2D moving
    mov_first5000 = squeeze(niftiread(fullfile(geneexpfolder, curgeneID))); % Gene data

    % Visual QC check
    figure
    imshowpair(movimg, mov_first5000); % Pre-registration overlay
    
    % Load pre-computed transformation
    load(strcat(check_dir, matfiles(i).name)); 
    movresult = movresult.movresult; % Contains transformation matrices
    
    % Apply spatial transformation
    regresult = imwarp(mov_first5000, movresult.Transformation,...
        'OutputView', imref2d(size(fiximg)));
    
    % Configure output directory
    cursavepath = fullfile(savepath, 'reg_gene_cluster', curID);
    if ~exist(cursavepath, 'dir')
        mkdir(cursavepath)
    end
    
    % Save transformed gene data
    for j = 1:size(regresult,3)
        niftiwrite(single(regresult(:,:,j)),...
            fullfile(cursavepath, sprintf('%s-%05d', curID, j)),...
            'Compressed', true);
    end
    
    clear regresult % Memory management
end
