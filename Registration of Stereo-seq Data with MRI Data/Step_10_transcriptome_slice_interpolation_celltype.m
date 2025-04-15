% Add SPM12 toolboxes for NIfTI processing
addpath /mnt/Data16Ta/Data/xiaojia/spm12
addpath /mnt/HotBackup/scratch/xiaojia/data/APP/spm12

%% === Convert cell type dot maps to density maps ===
layer_num = 6;
cell_type_num = 207;

% Load cortical layer mask
cortex_layer_name = '/mnt/SSD8T/xiaojia/sm_head_cortex_LIP2_head_T-all_to_MRI.nii.gz';
cortex_layer_data = spm_vol(cortex_layer_name);
cortex_layer_fdata = spm_read_vols(cortex_layer_data);

% Initialize cell array storage
cell_data = cell(1, cell_type_num);

% Load all cell type data
for i = 1:cell_type_num
    folder_path = '/mnt/SSD8T/xiaojia/cell_label_to_celltype_mask_iso10_to_MRI_merge/';
    file_name = strcat('head_LIP2_head_cell_type_', sprintf('%03d', i), '_sm_T-all_to_MRI.nii.gz');
    full_file_path = fullfile(folder_path, file_name);
    
    % Load and store cell type data
    small125_MRI_data = spm_vol(full_file_path);
    fdata2 = spm_read_vols(small125_MRI_data);
    cell_data{i} = fdata2;
end

% Save raw cell data
save('/mnt/SSD8T/xiaojia/cell_label_to_celltype_mask_iso10_to_MRI_merge/cell_data.mat', 'cell_data', '-v7.3');

%% === Density map processing ===
% Create 3D Gaussian filter for smoothing
sigma = 1;
kernel_size = 8;
gaussian_filter = fspecial3('gaussian', [kernel_size, kernel_size, 2], sigma);

% Initialize storage arrays
original_data_in_mask = cell(layer_num, cell_type_num);
density_maps = cell(layer_num, cell_type_num);

% Process each layer and cell type
for layer = 1:layer_num
    % Create layer-specific mask
    layer_mask = cortex_layer_fdata;
    layer_mask(layer_mask ~= layer) = 0;
    layer_mask(layer_mask == layer) = 1;
    
    for i = 1:cell_type_num
        % Apply layer mask to original data
        original_data_in_mask{layer, i} = cell_data{i} .* layer_mask;
        
        % Generate density maps through 3D convolution
        density_maps{layer, i} = convn(original_data_in_mask{layer, i}, gaussian_filter, 'same');
    end
end

% Sum density maps across layers
density_maps_sum = cell(1, cell_type_num);
for j = 1:cell_type_num
    temp_sum = density_maps{1, j};
    for i = 2:layer_num
        temp_sum = temp_sum + density_maps{i, j};
    end
    density_maps_sum{j} = temp_sum;
end

%% === Visualization ===
figure;
for i = 1:10
    % Original data visualization
    subplot(4, 5, i);
    slice_original = cell_data{i}(:, :, 70);
    imagesc(slice_original);
    colorbar;
    title(['Original Cell Type ', num2str(i)]);
    
    % Density map visualization
    subplot(4, 5, i + 10);
    slice_density = density_maps_sum{i}(:, :, 70);
    imagesc(slice_density);
    colorbar;
    title(['Fitted Cell Type ', num2str(i)]);
end

%% === Save density maps ===
for i = 1:cell_type_num
    output_name = strcat('/cell_label_to_celltype_mask_iso10_to_MRI_merge_density/',...
                        'head_LIP2_head_cell_type_', sprintf('%03d', i), 'sm_T-all_to_MRI.nii');
    cortex_layer_data.fname = output_name;
    spm_write_vol(cortex_layer_data, density_maps_sum{i});
end

%% === Spatial interpolation processing ===
% Load high-resolution cortex mask
cortex_layer_name = '/marmoset_gene_cortex/kj_ants_to_mri_2D/MBM_cortex_mask_25um_L_rotate_x6_z1.nii.gz';
cortex_mask_data = spm_vol(cortex_layer_name);
cortex_mask_fdata = spm_read_vols(cortex_mask_data);

% Prepare NaN mask for interpolation
cortex_mask_fdata_NaN = cortex_mask_fdata;
cortex_mask_fdata_NaN(cortex_mask_fdata_NaN ~= 0) = NaN;

% Load slice indices
slice125_num = load('/marmoset_gene_cortex/kj_ants_to_mri_2D/T402-T528_iso0025');

% Process each cell type
for i = 1:cell_type_num
    tic;
    % Load density map
    input_file = strcat('/cell_label_to_celltype_mask_iso10_to_MRI_merge_density/',...
                       'head_LIP2_head_cell_type_', sprintf('%03d', i), 'sm_T-all_to_MRI.nii');
    fdata2 = spm_read_vols(spm_vol(input_file));
    
    % Data orientation adjustment
    fdata2 = flip(permute(fdata2, [1, 3, 2]), 1);
    
    % Mask application
    masked_data = cortex_mask_fdata_NaN;
    masked_data(:, slice125_num, :) = fdata2;
    masked_data = masked_data .* cortex_mask_fdata;
    
    % Linear interpolation
    filled_data = fillmissing(masked_data, 'linear', 2, 'EndValues', 'nearest');
    
    % Save result
    output_name = strcat('/cell_label_to_celltype_mask_iso10_to_MRI_merge_density_chazhi/',...
                        'head_LIP2_head_cell_type_', sprintf('%03d', i), 'sm_T-all_to_MRI.nii');
    cortex_mask_data.fname = output_name;
    spm_write_vol(cortex_mask_data, filled_data);
    
    % Clear temporary data
    clear masked_data filled_data;
    toc;
end
