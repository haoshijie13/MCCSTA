%% Image Registration Quality Check Script
% Purpose: Verify alignment quality between registered 3D RNA data and original 2D slices

%% File Path Configuration
fixfolder = '/mnt/.../T402-528_resample01_bin5like'; % Registered 3D data path
movfolder = '/mnt/.../T402-528_resample01bin5like';   % Original 2D data path
check_dir = '/mnt/.../coronal_resample_bin5like/';   % Quality check output path

%% File List Initialization
fixfiles = dir(fullfile(fixfolder, '*.nii.gz'));  % Get registered files
movfiles = dir(fullfile(movfolder, '*.nii.gz'));   % Get original files

%% Main Processing Loop
for i = 1:numel(movfiles)
    %% File Loading
    curfixID = fixfiles(i).name;
    curID = curfixID(1:4);  % Extract slice ID (e.g., T402)
    curmovID = movfiles(i).name;
    
    % Read NIfTI files
    fiximg = squeeze(niftiread(fullfile(fixfolder, curfixID)));  % Registered 3D
    movimg = squeeze(niftiread(fullfile(movfolder, curmovID)));  % Original 2D

    %% Intensity Normalization
    % Slice-specific intensity clipping to optimize visualization
% ============ Intensity Adjustment Blocks ============

        if i == 5
            fiximg(fiximg >= 0.85) = 0.85;
            movimg(movimg >= 1800)  = 1800;
        end

        if i == 12
            fiximg(fiximg >= 0.75) = 0.75;
            movimg(movimg >= 1500) = 1500;
        end
 
        if i == 17
            fiximg(fiximg >= 0.83) = 0.83;
            movimg(movimg >= 2600) = 2600;
        end

        if i == 18
            fiximg(fiximg >= 0.8)  = 0.8;
            movimg(movimg >= 3300) = 3300;
        end

        if i == 19
            fiximg(fiximg >= 0.9)  = 0.9;
            movimg(movimg >= 4000) = 4000;
        end

        if i == 22
            fiximg(fiximg >= 0.4)  = 0.4;
            movimg(movimg >= 6000) = 6000;
        end

        if i == 28
            fiximg(fiximg >= 0.35) = 0.35;
            movimg(movimg >= 3050) = 3050;
        end

        if i == 24
            fiximg(fiximg >= 0.6)  = 0.6;
            movimg(movimg >= 5000) = 5000;
        end

        if i == 31
            fiximg(fiximg >= 0.9)  = 0.9;
            movimg(movimg >= 6500) = 6500;
        end

        if i == 44
            fiximg(fiximg >= 0.5)  = 0.5;
            movimg(movimg >= 6000) = 6000;
        end

        if i == 48
            fiximg(fiximg >= 0.8)   = 0.8;
            movimg(movimg >= 11000) = 11000;
        end

        if i == 49
            fiximg(fiximg >= 0.75) = 0.75;
            movimg(movimg >= 8500) = 8500;
        end

        if i == 52
            fiximg(fiximg >= 0.3) = 0;
            movimg(movimg >= 950) = 0;
 
        end

        if i == 56
            fiximg(fiximg >= 0.5)  = 0.5;
            movimg(movimg >= 6000) = 6000;
        end

        if i == 61
            fiximg(fiximg >= 0.82) = 0.82;
            movimg(movimg >= 3000) = 3000;
        end

        if i == 68
            fiximg(fiximg >= 0.8)   = 0.8;
            movimg(movimg >= 12000) = 12000;
        end
        
        if i == 86
            fiximg(fiximg >= 0.8)  = 0.8;
            movimg(movimg >= 9000) = 9000;
        end

        if i == 90
            fiximg(fiximg >= 0.5)  = 0.5;
            movimg(movimg >= 4000) = 4000;
        end

        if i == 101
            fiximg(fiximg >= 0.4) = 0.4;
            movimg(movimg >= 8000) = 0;
        end

        if i == 103
            fiximg(fiximg >= 0.75) = 0.75;
            movimg(movimg >= 5000) = 5000;
        end

        if i == 106
            fiximg(fiximg >= 0.4) = 0.4;
            movimg(movimg >= 3000) = 3000;
        end

        if i == 107
            fiximg(fiximg >= 0.4) = 0.4;
            movimg(movimg >= 3000) = 3000;
        end

        if i == 109
            fiximg(fiximg >= 0.5)  = 0.5;
            movimg(movimg >= 3500) = 3500;
        end

        if i == 110
            fiximg(fiximg >= 0.4) = 0.4;
            movimg(movimg >= 3000) = 3000;
        end

        if i == 111
            fiximg(fiximg >= 0.5)  = 0.5;
            movimg(movimg >= 2500) = 2500;
        end

        if i == 112
            fiximg(fiximg >= 0.85) = 0.85;
            movimg(movimg >= 4500) = 4500;
        end

        if i == 114
            fiximg(fiximg >= 0.3)   = 0.3;
            movimg(movimg >= 10000) = 10000;
        end
% ============ Visualization Block ============

obj = imshowpair(fiximg, fiximg);
figure
obj = imshowpair(movimg, movimg);

% ============ Final Registration ============
movresult = registerImages(movimg,fiximg);
save(strcat(check_dir , curID,'movresult.mat'),'movresult');

% ============ Result Export ============
saveas(obj, [check_dir , curID, 'fix.png'])        
obj = imshowpair(fiximg, movresult.RegisteredImage);
saveas(obj, [check_dir , curID, 'fixresult.png'])        
obj = imshowpair(fiximg, movresult.RegisteredImage,'montage');
saveas(obj, [check_dir ,curID, 'both.png'])
        

    %% Core Registration
    movresult = registerImages(movimg, fiximg);  % Execute image registration
    
    %% Result Saving
    % 1. Save transformation matrix
    save(fullfile(check_dir, [curID 'movresult.mat']), 'movresult');
    
    % 2. Generate verification images
    % Original vs Original (baseline)
    fig1 = imshowpair(fiximg, fiximg);
    saveas(fig1, fullfile(check_dir, [curID '_fix.png']));
    
    % Registered vs Original overlay
    fig2 = imshowpair(fiximg, movresult.RegisteredImage);
    saveas(fig2, fullfile(check_dir, [curID '_overlay.png']));
    
    % Side-by-side comparison
    fig3 = imshowpair(fiximg, movresult.RegisteredImage, 'montage');
    saveas(fig3, fullfile(check_dir, [curID '_comparison.png']));
end
