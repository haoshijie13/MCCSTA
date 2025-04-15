%% Image Registration Quality Control Pipeline
% This script performs visual verification of 3D-2D registration results
% between processed RNA expression maps and original spatial data

%% ===================== PATH CONFIGURATION ===========================
% Fixed images (registered 3D data)
fixfolder = '/marmoset_V2_20220627/coronal_2Dslice_nifti/large_rotate_reg_to_3D/totalRNA_brainmasked/T404-519';

% Moving images (original 2D slices)
movfolder = '/marmoset_V2_20220627/coronal_2Dslice_nifti/small/totalRNA_brainmasked/T404-519';   

% Output directory for QC results
check_dir = '/marmoset_V2_20220627_nifti/coronal/reg_check/Step2_5_register_headsurftorawsmall_check_coronal_check/';

%% ===================== FILE PREPARATION ==============================
% Get file lists using directory scanning
fixfiles = dir(fullfile(fixfolder, '*.nii.gz')); % Registered 3D files
movfiles = dir(fullfile(movfolder, '*.nii.gz'));  % Original 2D files

%% ===================== MAIN PROCESSING LOOP =========================
for i = 1:116  % Process all 116 slices
    
    %% ========== FILE IDENTIFICATION & LOADING ==========
    % Get current file IDs
    curfixID = fixfiles(i).name;
    curID = curfixID(1:4);  % Extract 4-character ID (e.g., T404)
    disp(['Processing: ' curID])
    curmovID = movfiles(i).name;
    
    % Construct full file paths
    fixfile = fullfile(fixfolder, curfixID);   % Registered 3D file
    movfile = fullfile(movfolder, curmovID);   % Original 2D file
    
    %% ========== IMAGE LOADING & PREPROCESSING ==========
    % Read NIfTI files and remove singleton dimensions
    movimg = squeeze(niftiread(movfile));  % Original RNA data
    fiximg = squeeze(niftiread(fixfile));  % Registered RNA data
    
    %% ========== INTENSITY NORMALIZATION ================
    % Case-specific intensity clipping to optimize visualization
    % Note: Thresholds determined empirically for each slice
    switch i
        case 5  % Slice 5 adjustments
            fiximg(fiximg >= 0.85) = 0.85;
            movimg(movimg >= 1800) = 1800;
            
        case 12  % Slice 12 adjustments
            fiximg(fiximg >= 0.75) = 0.75;
            movimg(movimg >= 1500) = 1500;
            
        case 17             
            % Adjust slice 17 intensity ranges
            fiximg(fiximg >= 0.83) = 0.83;  % Clip high values in registered data
            movimg(movimg >= 2600) = 2600;   % Clip original data intensities

        case 18             
            % Optimize contrast for slice 18
            fiximg(fiximg >= 0.8) = 0.8;    % Upper threshold for registration result
            movimg(movimg >= 3300) = 3300;  % Corresponding raw data adjustment

        case 19            
            % Handle saturation in slice 19
            fiximg(fiximg >= 0.9) = 0.9;    % Prevent overexposure in 3D data
            movimg(movimg >= 4000) = 4000;  % Match 2D data dynamic range

        case 22             
            % Balance low-contrast features in slice 22
            fiximg(fiximg >= 0.4) = 0.4;    % Reduce background noise
            movimg(movimg >= 6000) = 6000;  % Preserve structural details

        case 24            
            % Enhance mid-range features for slice 24
            fiximg(fiximg >= 0.6) = 0.6;    % Optimize mid-tone visibility
            movimg(movimg >= 5000) = 5000;  % Match intensity profiles

        case 28             
            % Improve edge visibility in slice 28
            fiximg(fiximg >= 0.35) = 0.35;  % Sharpen boundary definition
            movimg(movimg >= 3050) = 3050;  % Maintain tissue contrast

        case 31
            % Handle extreme values in slice 31
            fiximg(fiximg >= 0.9) = 0.9;    % Control registration artifacts
            movimg(movimg >= 6500) = 6500;  % Normalize intensity outliers

        case 44              
            % Adjust for uneven illumination in slice 44
            fiximg(fiximg >= 0.5) = 0.5;    % Compensate light falloff
            movimg(movimg >= 6000) = 6000;  % Balance intensity distribution

        case 49              
            % Manage high dynamic range in slice 49
            fiximg(fiximg >= 0.75) = 0.75;  % Preserve subtle intensity variations  
            movimg(movimg >= 8500) = 8500;  % Prevent signal clipping

        case 52             
            % Special handling for low-signal slice 52
            fiximg(fiximg >= 0.3) = 0;      % Remove residual background
            movimg(movimg >= 950) = 0;      % Eliminate noise floor

        case 56             
            % Enhance cortical layer contrast in slice 56
            fiximg(fiximg >= 0.5) = 0.5;    % Highlight laminar organization
            movimg(movimg >= 6000) = 6000;  % Match histological features

        case 61            
            % Optimize hippocampus visualization in slice 61
            fiximg(fiximg >= 0.85) = 0.85;  % Improve subfield contrast
            movimg(movimg >= 3000) = 3000;  % Maintain anatomical fidelity

        case 66             
            % Adjust for partial volume effects in slice 66
            fiximg(fiximg >= 0.2) = 0.2;    % Reduce edge blurring
            movimg(movimg >= 1800) = 1800;  % Enhance structural clarity

        case 68            
            % Handle hyperintense regions in slice 68
            fiximg(fiximg >= 0.8) = 0.8;    % Control vessel signal
            movimg(movimg >= 12000) = 12000; % Prevent blooming artifacts

        case 75             
            % Improve white matter contrast in slice 75
            fiximg(fiximg >= 0.45) = 0.45;  % Differentiate WM/GM
            movimg(movimg >= 6000) = 6000;  % Maintain tract visibility

        case 90             
            % Standardize intensity for comparative analysis
            fiximg(fiximg >= 0.5) = 0.5;    % Normalize across specimens
            movimg(movimg >= 4000) = 4000;  % Cross-sample calibration
        
        case 114  % Final slice adjustments
            fiximg(fiximg >= 0.3) = 0.3;
            movimg(movimg >= 10000) = 10000;
    end
    
    %% ========== IMAGE REGISTRATION =====================
    % Core registration using MATLAB's registerImages function
    movresult = registerImages(movimg, fiximg);
    
    % Save transformation matrix for potential reuse
    save(fullfile(check_dir, [curID '_movresult.mat']), 'movresult');
    
    %% ========== VISUALIZATION & OUTPUT =================
    % 1. Fixed Image Visualization
    fig1 = imshowpair(fiximg, fiximg);
    saveas(fig1, fullfile(check_dir, [curID '_fix.png']));
    
    % 2. Registration Result Overlay
    fig2 = imshowpair(fiximg, movresult.RegisteredImage);
    saveas(fig2, fullfile(check_dir, [curID '_overlay.png']));
    
    % 3. Side-by-Side Comparison
    fig3 = imshowpair(fiximg, movresult.RegisteredImage, 'montage');
    saveas(fig3, fullfile(check_dir, [curID '_comparison.png']));
    
end
