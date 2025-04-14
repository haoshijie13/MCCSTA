%% Segment gene/celltypes and networks using PLS or correlation

% Load data matrices
savedir = '/data/gene/result/marmoset/all_map_fit/fit_by_gene_segment/';
data1 = readmatrix('/data/gene/result/marmoset/hz_result/Segment_network_matrix.csv');
data1 = data1(:,2:end); % Remove index column

data2 = readmatrix('/data/gene/result/marmoset/hz_result/Segment_gene_matrix.csv');
data2 = data2(:,2:end); % Remove index column

num_components = 20; % Number of PLS components

%% Perform PLS regression for each network
for net = 1:size(data1,2)
    % Set up PLS inputs
    x = data2;
    y = data1(:,net);
    
    % Run PLS regression
    [XL, YL, XS, YS, BETA, PCTVAR, MSE, stats] = plsregress(x, y, num_components);
    
    % Calculate variable contributions
    variable_contribution = sum(stats.W.^2, 2);
    
    % Generate model predictions
    y_fit = [ones(size(x,1),1) x] * BETA;
    
    % Calculate model performance
    [corr_p, corr_q] = corr(y, y_fit);
    
    % Save regression coefficients
    save(fullfile(savedir, ['marmoset_net',num2str(net,'%02d'),...
        '_PLS' num2str(num_components) '_y_fit_BETA.mat']), 'y_fit', 'BETA');
end

%% Export PLS results to CSV
data_dir = dir('/data/gene/result/marmoset/all_map_fit/fit_by_gene_segment/net_PLS/marmoset_net*');
b = [];

% Aggregate beta coefficients
for i = 1:15
    data = load(fullfile(data_dir(i).folder, data_dir(i).name));
    b = [b, data.BETA];
end

% Prepare column headers
network_table = readtable('/data/gene/result/marmoset/hz_result/Segment_network_matrix.csv');
gene_table = readtable('/data/gene/result/marmoset/hz_result/Segment_gene_matrix.csv');

col_names = [{'BETA'}, gene_table.Properties.VariableNames(2:end)];

% Create and save results table
b_table = array2table(b, 'VariableNames', network_table.Properties.VariableNames(2:end));
b_table.Properties.RowNames = col_names;
writetable(b_table, fullfile(savedir, 'marmoset_net_PLS20_y_fit_BETA.csv'));

%% Gene contribution analysis
% Load PLS results
data_dir = dir(fullfile(savedir, 'net_PLS/*.mat'));
data_PLS_all = [];

% Aggregate all network coefficients
for net = 1:15
    data = load(fullfile(data_dir(net).folder, data_dir(net).name));
    data_PLS_all = [data_PLS_all, data.BETA];
end

% Set analysis parameters
x = data2;
top_gene_counts = 1:5000; % Full range of genes to test
num_networks = 15;
corr_thr = 0.95;

% Initialize result matrices
results = zeros(num_networks, length(top_gene_counts));
y_fit_all = zeros(num_networks, length(top_gene_counts), size(data1, 1));
top_counts_at_thr = zeros(num_networks, 1);

% Analyze gene contributions per network
for net = 1:num_networks
    BETA = data_PLS_all(:, net);
    gene_weights = BETA(2:end); % Exclude intercept
    
    % Test different gene subset sizes
    for top_idx = 1:length(top_gene_counts)
        top_count = top_gene_counts(top_idx);
        
        % Select top contributing genes
        [~, top_indices] = maxk(abs(gene_weights), top_count);
        x_selected = x(:, top_indices);
        BETA_selected = [BETA(1); gene_weights(top_indices)];
        
        % Calculate predictions
        y_fit = [ones(size(x_selected, 1), 1), x_selected] * BETA_selected;
        y_fit_all(net, top_idx, :) = y_fit;
        
        % Store correlation results
        results(net, top_idx) = corr(y_fit, data1(:, net));
        
        % Check correlation threshold
        if results(net, top_idx) >= corr_thr
            top_counts_at_thr(net) = top_count;
            break;
        end
    end
end

% Save analysis results
save(fullfile(savedir, 'net_gene_PLS20_gene_num_and_fitting_results_thr095.mat'),...
    'results', 'top_counts_at_thr', '-v7.3');
writematrix(top_gene_counts, fullfile(savedir, 'top_gene_counts_thr095.csv'));
writematrix(top_counts_at_thr, fullfile(savedir, 'top_counts_at_thr_thr095.csv'));
