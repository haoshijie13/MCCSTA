% Add SPM12 toolboxes for NIfTI file processing
addpath /mnt/Data18Tb/Data/xiaojia/spm12;
addpath /zfs/scratch/xiaojia/data/APP/spm12;

% Load cortex mask data (25um resolution, left hemisphere)
cortex_mask_name = '/marmoset_gene_cortex/kj_ants_to_mri_2D/MBM_cortex_mask_25um_L_rotate_x6_z1.nii.gz';
cortex_mask_data = spm_vol(cortex_mask_name);
cortex_mask_fdata = spm_read_vols(cortex_mask_data);

% Prepare NaN mask for missing value handling
cortex_mask_fdata_NaN = cortex_mask_fdata;
cortex_mask_fdata_NaN(cortex_mask_fdata_NaN ~= 0) = NaN; % Convert mask areas to NaN

% Load slice indices for coronal plane processing
slice125_num = load('/marmoset_gene_cortex/kj_ants_to_mri_2D/T402-T528_iso0025');

% Main processing loop for gene data slices 4451-4500
for i = 4451:4500
    % Start overall timer
    totalTimer = tic;
    
    % Check file existence before processing
    folder_path = '/marmoset_coronal_gene/sevent5000_merge_3D_merge_chazhi/';
    file_name = strcat('gene_', sprintf('%04d', i), '_in_MRI_chazhi.nii');
    full_file_path = fullfile(folder_path, file_name);
    
    if exist(full_file_path, 'file')
        disp(['File ', file_name, ' already exists. Skipping...']);
        continue;
    end

    % Load gene data and preprocess
    small125_MRI_data_name = strcat('/coronal_gene_to_MRI/apply_result_gene/sevent5000_merge_3D_merge/gene_', sprintf('%04d', i), '_in_MRI.nii.gz');
    small125_MRI_data = spm_vol(small125_MRI_data_name);
    fdata2 = spm_read_vols(small125_MRI_data);
    
    % Permute and flip data for coronal plane alignment
    fdata2 = permute(fdata2, [1, 3, 2]);  % Adjust axis orientation
    fdata2 = flip(fdata2, 1);            % Flip left-right for hemisphere alignment

    % Apply gene data to cortex mask
    cortex_mask_fdata_NaN_eachgene = cortex_mask_fdata_NaN;
    cortex_mask_fdata_NaN_eachgene(:, slice125_num, :) = fdata2;
    cortex_mask_fdata_NaN_eachgene = cortex_mask_fdata_NaN_eachgene .* cortex_mask_fdata; % Mask application

    % Linear interpolation for missing values (NaN)
    fillmiss_cortex_mask_fdata = fillmissing(cortex_mask_fdata_NaN_eachgene, 'linear', 2, 'EndValues', 'nearest');
    
    % Save processed data
    cortex_mask_data_eachgene = cortex_mask_data;
    cortex_mask_data_eachgene.fname = full_file_path;
    spm_write_vol(cortex_mask_data_eachgene, fillmiss_cortex_mask_fdata);
    
    % Clear temporary variables
    clear cortex_mask_fdata_NaN_eachgene fillmiss_cortex_mask_fdata;
    
    % Display timing information
    elapsedTime = toc(totalTimer);
    disp(['Processing completed in: ', num2str(elapsedTime), ' seconds']);
end

% Clean workspace after processing
clear;
