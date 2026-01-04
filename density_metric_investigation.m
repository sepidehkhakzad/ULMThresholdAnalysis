%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Ultrasound Super-Resolution Imaging: Error Analysis Tool for different densities.
%
% Author: Sepideh K. Gharamaleki
% Lab: IMPACT LAB, Dr. Hassan Rivaz
% Date: 2026/01/03
% Version: 1.0
%
% Description:
% This script performs sensitivity analysis for ultrasound super-resolution
% imaging by introducing controlled false positives and false negatives
% into ground truth localization data. It evaluates the impact on image
% quality metrics (SSIM, Dice coefficient, PSNR) for different density regions.
%
% Citation:
% If you use this code in your research, please cite:
% [Sepideh K. Gharamaleki et al., "Evaluating Detection Thresholds: The Impact of False Positives and Negatives on Super-Resolution Ultrasound Localization Microscopy," SPIE Medical Imaging, 2026]
%
% Requirements:
% - MATLAB
% - Image Processing Toolbox
% - Parallel Computing Toolbox (for parfor loops)
% - Input files: gt_scat_inside2.txt, metadata_simu2.mat Downloadable at
%       the link provided in the github repo. 
%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

close all
clear

% Load ground truth data and metadata
ground_truth = load('gt_scat_inside2.txt');
load('metadata_simu2.mat');
x_grid = repmat((PxSet.mapX), 512, 1);
z_grid = repmat((PxSet.mapZ)', 1, 512);

cf = SimSet.centre_frequency;
% Define the metrics
lambdaa = 1540/cf;
threshold = lambdaa/2;   % threshold for matching
threshold_regions = 3000; %3000 for Simu2 and 800 for Simu1 experimentally computed
grid_size = [512, 512];  % grid size is set based on the resolution
gt_SRDensity = compute_gt_density_map(ground_truth, PxSet);

load("kde_density_simu2.mat");
dense_regions = classifyRegions(density, threshold_regions, "dense");
sparse_regions = classifyRegions(density, threshold_regions, "sparse");
% Show dense vs sparse regions, KDE and SR density map
figure;
subplot(2, 2, 1);
imagesc(gt_SRDensity'); colormap('hot'); title('SR Density Map');
subplot(2, 2, 2);
imagesc(density); colormap('hot'); title('KDE Density Map');
subplot(2, 2, 3);
imagesc(sparse_regions); title('Sparse Regions');
subplot(2, 2, 4);
imagesc(dense_regions); title('Dense Regions');



% Define the ranges for num_fp_per_frame and fn_rate
num_fp_range = 0:10:100; 
fn_rate_range = 0:0.1:0.8; 
fp_rate_range = 0:0.1:0.8;
num_fn_range = 0:10:100;

% Initialize matrices to store the metrics
ssim_matrix_number = zeros(length(num_fp_range), length(num_fn_range), 4);
dice_matrix_number = zeros(length(num_fp_range), length(num_fn_range), 4);
psnr_matrix_number = zeros(length(num_fp_range), length(num_fn_range), 4);

ssim_matrix_rate = zeros(length(fp_rate_range), length(fn_rate_range), 4);
dice_matrix_rate = zeros(length(fp_rate_range), length(fn_rate_range), 4);
psnr_matrix_rate = zeros(length(fp_rate_range), length(fn_rate_range), 4);

dense_regions = dense_regions';
sparse_regions = sparse_regions';

ssim_tmp = cell(1,8);
dice_tmp = cell(1,8);
psnr_tmp = cell(1,8);

% Iterate over the ranges and compute metrics
% This will take several hours.
parfor mode = 1:8
    disp(sprintf("mode = %d", mode))

    isNumber = mode <= 4; % whether the experiment is based on the number of MBs or percentage (i.e. rate)
    localIdx = mod(mode-1,4) + 1; % variable to determine whether the experiment is for dense or sparse regions.

    if localIdx <= 2
        region = 'dense';
        mask = dense_regions;
    else
        region = 'sparse';
        mask = sparse_regions;
    end

    isPerFrame = mod(localIdx,2) == 0; % Variable to determine whether the changes are per frame or for the whole dataset

    if isNumber
        ssim_local = zeros(length(num_fp_range), length(num_fn_range));
        dice_local = zeros(length(num_fp_range), length(num_fn_range));
        psnr_local = zeros(length(num_fp_range), length(num_fn_range));
    else
        ssim_local = zeros(length(fp_rate_range), length(fn_rate_range));
        dice_local = zeros(length(fp_rate_range), length(fn_rate_range));
        psnr_local = zeros(length(fp_rate_range), length(fn_rate_range));
    end

    if isNumber
        for i = 1:length(num_fp_range)
            for j = 1:length(num_fn_range)
                num_fp = num_fp_range(i);
                num_fn = num_fn_range(j);

                if isPerFrame
                    modified_data = modifyDataPerFrameNumberSpecificRegion( ...
                        ground_truth, x_grid, z_grid, ...
                        num_fn, num_fp, threshold, mask, region);
                else
                    modified_data = modifyDataTotalNumberSpecificRegion( ...
                        ground_truth, x_grid, z_grid, ...
                        num_fn, num_fp, threshold, mask, region);
                end

                modified_SRDensity = compute_gt_density_map(modified_data, PxSet);
                mask_resize = imresize(mask, [size(mask,1)*10 size(mask,2)*10]);

                ssim_local(i,j) = calculate_ssim_mask(gt_SRDensity, modified_SRDensity, mask_resize);
                dice_local(i,j) = calculate_dice_mask(gt_SRDensity, modified_SRDensity, mask_resize);
                psnr_local(i,j) = calculate_psnr_mask(gt_SRDensity, modified_SRDensity, mask_resize);
            end
        end
    else
        for i = 1:length(fp_rate_range)
            for j = 1:length(fn_rate_range)
                fp_rate = fp_rate_range(i);
                fn_rate = fn_rate_range(j);

                if isPerFrame
                    modified_data = modifyDataPerFrameRateSpecificRegion( ...
                        ground_truth, x_grid, z_grid, ...
                        fn_rate, fp_rate, threshold, mask, region);
                else
                    modified_data = modifyDataTotalRateSpecificRegion( ...
                        ground_truth, x_grid, z_grid, ...
                        fn_rate, fp_rate, threshold, mask, region);
                end

                modified_SRDensity = compute_gt_density_map(modified_data, PxSet);
                mask_resize = imresize(mask, [size(mask,1)*10 size(mask,2)*10]);

                ssim_local(i,j) = calculate_ssim_mask(gt_SRDensity, modified_SRDensity, mask_resize);
                dice_local(i,j) = calculate_dice_mask(gt_SRDensity, modified_SRDensity, mask_resize);
                psnr_local(i,j) = calculate_psnr_mask(gt_SRDensity, modified_SRDensity, mask_resize);
            end
        end
    end

    ssim_tmp{mode} = ssim_local;
    dice_tmp{mode} = dice_local;
    psnr_tmp{mode} = psnr_local;
end

% assignments to final matrices
for mode = 1:8
    localIdx = mod(mode-1,4) + 1;

    if mode <= 4
        ssim_matrix_number(:,:,localIdx) = ssim_tmp{mode};
        dice_matrix_number(:,:,localIdx) = dice_tmp{mode};
        psnr_matrix_number(:,:,localIdx) = psnr_tmp{mode};
    else
        ssim_matrix_rate(:,:,localIdx) = ssim_tmp{mode};
        dice_matrix_rate(:,:,localIdx) = dice_tmp{mode};
        psnr_matrix_rate(:,:,localIdx) = psnr_tmp{mode};
    end
end

% The following visualizes the computed matrix
% There are 8 experiments, 3 metrics for each experiment,
% yeilding 24 plots.

[fp_rate_range_grid, fn_rate_grid] = meshgrid(fp_rate_range, fn_rate_range);
[num_fp_range_grid, num_fn_range_grid] = meshgrid(num_fp_range, num_fn_range);

% Each entry = one experimental condition
experiments = struct( ...
    'label',  {'N-T-D','N-PF-D','N-T-S','N-PF-S','R-T-D','R-PF-D','R-T-S','R-PF-S'}, ...
    'type',   {'Number','Number','Number','Number','Rate','Rate','Rate','Rate'}, ...
    'index',  {1,2,3,4,1,2,3,4}, ...
    'mode',   {'Total','Per Frame','Total','Per Frame', ...
               'Total','Per Frame','Total','Per Frame'}, ...
    'region', {'Dense','Dense','Sparse','Sparse', ...
               'Dense','Dense','Sparse','Sparse'} ...
);

metrics = {'SSIM', 'Dice', 'PSNR'};

for e = 1:numel(experiments)
    exp = experiments(e);

    % Axis + data source chosen once per experiment
    if strcmp(exp.type,'Number')
        X = num_fp_range_grid;
        Y = num_fn_range_grid;
        Data = {ssim_matrix_number, dice_matrix_number,psnr_matrix_number};
    else
        X = fp_rate_range_grid;
        Y = fn_rate_grid;
        Data = {ssim_matrix_rate, dice_matrix_rate, psnr_matrix_rate};
    end

    for m = 1:numel(metrics)

        % Metric selection
        switch metrics{m}
            case 'SSIM'
                Z = Data{1}(:,:,exp.index);
            case 'Dice'
                Z = Data{2}(:,:,exp.index);
            case 'PSNR'
                Z = Data{3}(:,:,exp.index);
        end

        % Plot
        plot_with_label( ...
            Z, X, Y, ...
            metrics{m}, ...
            exp.mode, ...
            exp.type, ...
            exp.region ...
        );
    end
end


function plot_with_label(thing_to_plot, x, y, what_is_being_plotted, total_or_perframe, number_or_rate, sparse_or_dense)
    
    figure;
    surf(x, y, thing_to_plot);
    xlabel("False Positive " + total_or_perframe + " " + number_or_rate );
    ylabel("False Negative " + total_or_perframe + " " + number_or_rate );
    zlabel(what_is_being_plotted);
    title(what_is_being_plotted+ " vs. False Negative and False Positive " + total_or_perframe + " " + number_or_rate+ " " + sparse_or_dense);
    colorbar;
end

function SR_density_resized = compute_gt_density_map(gt, PxSet)
    % Initialize grids and parameters
    x_grid = repmat((PxSet.mapX), 512, 1);
    z_grid = repmat((PxSet.mapZ)', 1, 512);
    x_steps = PxSet.mapX;
    z_steps = PxSet.mapZ;
    dx = PxSet.dx;
    dz = PxSet.dz;
    sigma_x = dx;
    sigma_z = dz;

    % Initialize the density matrix
    SR_density = zeros(512, 512);

    % Load ground truth scatterer positions
    indices_total = gt;

    % Normalize the scatterer positions with respect to the grid
    a = (indices_total(:, 3) - min(PxSet.mapX)) / dx;
    b = (indices_total(:, 5) - min(PxSet.mapZ)) / dz;
    indices_total = [a, b];

    % Compute the SR density map
    for jj = 1:size(indices_total, 1)
        current_center = indices_total(jj, :);

        if current_center(1) < 1 || current_center(2) < 1
            continue;
        end

        c_z = z_steps(floor(current_center(1))) + (current_center(1) - floor(current_center(1))) * 3.9139e-5;
        c_x = x_steps(floor(current_center(2))) + (current_center(2) - floor(current_center(2))) * 3.9139e-5;

        SR_density = SR_density + exp(-(((x_grid - c_x) / sigma_x).^2 + ((z_grid - c_z) / sigma_z).^2) / 2);
    end
    SR_density_resized = imresize(SR_density, [size(SR_density, 1) * 10, size(SR_density, 2) * 10]);
    
end

function data_with_fp = addFalsePositivesPerFrameRateSpecificRegion(data, dataforsize, x_grid, z_grid, fp_percentage, distance_threshold, region_mask, region_type)
    % Adds false positives either in dense or sparse regions based on region_type.
    % region_mask: Binary mask where 1 indicates a region of interest (dense/sparse), 0 otherwise.
    % region_type: Either 'dense' or 'sparse'
    
    data_with_fp = data;
    frames = unique(data(:, 1)); % Unique frame numbers
    
    
    for i = 1:length(frames)
        frame_idx = frames(i);
        frame_data = dataforsize(dataforsize(:, 1) == frame_idx, :);
        
        % Extract the true positive coordinates (x, z) for the current frame
        true_x = frame_data(:, 3); % Assuming x-coordinates are in the 3rd column
        true_z = frame_data(:, 5); % Assuming z-coordinates are in the 5th column
        
        % Calculate the number of false positives to add for this frame
        num_frame_points = size(frame_data, 1);
        num_fp_for_frame = round(num_frame_points * fp_percentage);
        
        for j = 1:num_fp_for_frame
            false_positive_added = false;
            while ~false_positive_added
                % Choose random position in the grid
                rand_x = randi(size(x_grid, 2)); % Random column index
                rand_z = randi(size(z_grid, 1)); % Random row index
                
                % Get the x and z coordinates of the randomly selected grid point
                fp_x = x_grid(rand_z, rand_x);
                fp_z = z_grid(rand_z, rand_x);
                
                % Check if the random point is in the desired region (dense/sparse)
                if (region_type == "dense" && region_mask(rand_z, rand_x) == 1) || ...
                   (region_type == "sparse" && region_mask(rand_z, rand_x) == 1)
                    % Compute the distance to all true positives in this frame
                    distances = sqrt((true_x - fp_x).^2 + (true_z - fp_z).^2);

                    % If the minimum distance is greater than the threshold, accept this point
                    if all(distances > distance_threshold)
                        % Create a random false positive
                        false_positive = [frame_idx, max(frame_data(:, 2)) + j, fp_x, 0, fp_z, 0];

                        % Add the false positive to the data
                        data_with_fp = [data_with_fp; false_positive];
                        false_positive_added = true;
                    end
                end
            end
        end
    end
end

function data_with_fn = addFalseNegativesPerFrameRateSpecificRegion(data, fn_rate, region_mask, region_type, x_grid, z_grid)
    % Adds false negatives either in dense or sparse regions based on region_type.
    % region_mask: Binary mask where 1 indicates a region of interest (dense/sparse), 0 otherwise.
    % region_type: Either 'dense' or 'sparse'
    
    data_with_fn = [];
    frames = unique(data(:, 1));  % Unique frame numbers
    prev_size = 0;
    for i = 1:length(frames)
        frame_idx = frames(i);
        frame_data = data(data(:, 1) == frame_idx, :);
        
        % Extract the true positive coordinates (x, z) for the current frame
        true_x = frame_data(:, 3); % Assuming x-coordinates are in the 3rd column
        true_z = frame_data(:, 5); % Assuming z-coordinates are in the 5th column
        
        % Calculate the number of false negatives to introduce
        
        % Check if the points are in the desired region
        indices_to_remove = [];
        for j = 1:size(frame_data, 1)
            % Convert true positive (x, z) to grid coordinates
            [~, rand_x] = min(abs(x_grid(1, :) - true_x(j)));
            [~, rand_z] = min(abs(z_grid(:, 1) - true_z(j)));
            
            % Check if the point is in the desired region
            if (region_type == "dense" && region_mask(rand_z, rand_x) == 1) || ...
               (region_type == "sparse" && region_mask(rand_z, rand_x) == 1)
                indices_to_remove = [indices_to_remove; j ]; % Mark this index for removal
            end
        end
        num_to_remove = round(length(indices_to_remove) * fn_rate);
        if ~isempty(indices_to_remove)
                num_to_remove = min(length(indices_to_remove), num_to_remove);
                remove_idx = indices_to_remove(randperm(length(indices_to_remove), num_to_remove));
                frame_data(remove_idx,:) = [];
                data_with_fn = [data_with_fn; frame_data];
        else
            disp("not enough data in the region to remove")
        end
    end
end

function data_with_fp = addFalsePositivesTotalRateSpecificRegion(data, x_grid, z_grid, fp_percentage, distance_threshold, region_mask, region_type)
    % Adds false positives based on the total rate for dense or sparse regions.
    % region_mask: Binary mask where 1 indicates a region of interest (dense/sparse), 0 otherwise.
    % region_type: Either 'dense' or 'sparse'
    
    num_original_points = size(data, 1);
    total_num_fp = round(num_original_points * fp_percentage);
    
    data_with_fp = data;
    frames = unique(data(:, 1)); % Unique frame numbers
    
    for i = 1:total_num_fp
        false_positive_added = false;
        while ~false_positive_added
            % Choose a random frame
            rand_frame = frames(randi(length(frames)));
            
            % Get data for the chosen frame
            frame_data = data(data(:, 1) == rand_frame, :);
            
            % Extract the true positive coordinates (x, z) for the current frame
            true_x = frame_data(:, 3); % Assuming x-coordinates are in the 3rd column
            true_z = frame_data(:, 5); % Assuming z-coordinates are in the 5th column
            
            % Choose random position in the grid
            rand_x = randi(size(x_grid, 2)); % Random column index
            rand_z = randi(size(z_grid, 1)); % Random row index
            
            % Get the x and z coordinates of the randomly selected grid point
            fp_x = x_grid(rand_z, rand_x);
            fp_z = z_grid(rand_z, rand_x);
            
            % Check if the random point is in the desired region (dense/sparse)
            if (region_type == "dense" && region_mask(rand_z, rand_x) == 1) || ...
               (region_type == "sparse" && region_mask(rand_z, rand_x) == 1)
                % Compute the distance to all true positives in this frame
                distances = sqrt((true_x - fp_x).^2 + (true_z - fp_z).^2);
                
                % If the minimum distance is greater than the threshold, accept this point
                if all(distances > distance_threshold)
                    % Create a random false positive
                    false_positive = [rand_frame, max(frame_data(:, 2)) + 1, fp_x, 0, fp_z, 0];
                    
                    % Add the false positive to the data
                    data_with_fp = [data_with_fp; false_positive];
                    false_positive_added = true;
                end
            end
        end
    end
end

function data_with_fp = addFalsePositivesTotalNumberSpecificRegion(data, x_grid, z_grid, total_num_fp, distance_threshold, region_mask, region_type)
    % Adds a total number of false positives for dense or sparse regions.
    % region_mask: Binary mask where 1 indicates a region of interest (dense/sparse), 0 otherwise.
    % region_type: Either 'dense' or 'sparse'
    
    data_with_fp = data;
    frames = unique(data(:, 1)); % Unique frame numbers
    
    for i = 1:total_num_fp
        false_positive_added = false;
        while ~false_positive_added
            % Choose a random frame
            rand_frame = frames(randi(length(frames)));
            
            % Get data for the chosen frame
            frame_data = data(data(:, 1) == rand_frame, :);
            
            % Extract the true positive coordinates (x, z) for the current frame
            true_x = frame_data(:, 3); % Assuming x-coordinates are in the 3rd column
            true_z = frame_data(:, 5); % Assuming z-coordinates are in the 5th column
            
            % Choose random position in the grid
            rand_x = randi(size(x_grid, 2)); % Random column index
            rand_z = randi(size(z_grid, 1)); % Random row index
            
            % Get the x and z coordinates of the randomly selected grid point
            fp_x = x_grid(rand_z, rand_x);
            fp_z = z_grid(rand_z, rand_x);
            
            % Check if the random point is in the desired region (dense/sparse)
            if (region_type == "dense" && region_mask(rand_z, rand_x) == 1) || ...
               (region_type == "sparse" && region_mask(rand_z, rand_x) == 1)
                % Compute the distance to all true positives in this frame
                distances = sqrt((true_x - fp_x).^2 + (true_z - fp_z).^2);
                
                % If the minimum distance is greater than the threshold, accept this point
                if all(distances > distance_threshold)
                    % Create a random false positive
                    false_positive = [rand_frame, max(frame_data(:, 2)) + 1, fp_x, 0, fp_z, 0];
                    
                    % Add the false positive to the data
                    data_with_fp = [data_with_fp; false_positive];
                    false_positive_added = true;
                end
            end
        end
    end
end

function data_with_fp = addFalsePositivesPerFrameNumberSpecificRegion(data, x_grid, z_grid, num_fp_per_frame, distance_threshold, region_mask, region_type)
    % Adds a number of false positives per frame for dense or sparse regions.
    % region_mask: Binary mask where 1 indicates a region of interest (dense/sparse), 0 otherwise.
    % region_type: Either 'dense' or 'sparse'
    
    data_with_fp = data;
    frames = unique(data(:, 1));  % Unique frame numbers
    
    for i = 1:length(frames)
        frame_idx = frames(i);
        frame_data = data(data(:, 1) == frame_idx, :);
        
        % Extract the true positive coordinates (x, z) for the current frame
        true_x = frame_data(:, 3); % Assuming x-coordinates are in the 3rd column
        true_z = frame_data(:, 5); % Assuming z-coordinates are in the 5th column
        
        for j = 1:num_fp_per_frame
            false_positive_added = false;
            while ~false_positive_added
                % Choose random position in the grid
                rand_x = randi(size(x_grid, 2)); % Random column index
                rand_z = randi(size(z_grid, 1)); % Random row index
                
                % Get the x and z coordinates of the randomly selected grid point
                fp_x = x_grid(rand_z, rand_x);
                fp_z = z_grid(rand_z, rand_x);
                
                % Check if the random point is in the desired region (dense/sparse)
                if (region_type == "dense" && region_mask(rand_z, rand_x) == 1) || ...
                   (region_type == "sparse" && region_mask(rand_z, rand_x) == 1)
                    % Compute the distance to all true positives in this frame
                    distances = sqrt((true_x - fp_x).^2 + (true_z - fp_z).^2);
                    
                    % If the minimum distance is greater than the threshold, accept this point
                    if all(distances > distance_threshold)
                        % Create a random false positive
                        false_positive = [frame_idx, max(frame_data(:, 2)) + j, fp_x, 0, fp_z, 0];
                        
                        % Add the false positive to the data
                        data_with_fp = [data_with_fp; false_positive];
                        false_positive_added = true;
                    end
                end
            end
        end
    end
    end

function data_with_fn = addFalseNegativesTotalRateSpecificRegion(data, fn_rate, region_mask, region_type, x_grid, z_grid)
    % Adds false negatives based on the total rate for dense or sparse regions.
    % region_mask: Binary mask where 1 indicates a region of interest (dense/sparse), 0 otherwise.
    % region_type: Either 'dense' or 'sparse'
    
    data_with_fn = data;
    frames = unique(data(:, 1));  % Unique frame numbers
    total_points = size(data, 1);
    
    % Total number of false negatives to remove
    
    points_removed = 0; % Track how many points have been removed
    indices_to_remove = [];
    prev_size = 0;
    if fn_rate ~=0
        for i = 1:length(frames)
            frame_idx = frames(i);
            frame_data = data(data(:, 1) == frame_idx, :);
            
            % Extract the true positive coordinates (x, z) for the current frame
            true_x = frame_data(:, 3); % Assuming x-coordinates are in the 3rd column
            true_z = frame_data(:, 5); % Assuming z-coordinates are in the 5th column
            prev_size = prev_size + size(data(data(:, 1) == frame_idx-1, :),1);
            % Check if the points are in the desired region
            
            for j = 1:size(frame_data, 1)
                % Convert true positive (x, z) to grid coordinates
                [~, rand_x] = min(abs(x_grid(1, :) - true_x(j)));
                [~, rand_z] = min(abs(z_grid(:, 1) - true_z(j)));
                
                % Check if the point is in the desired region (dense or sparse)
                if (region_type == "dense" && region_mask(rand_z, rand_x) == 1) || ...
                   (region_type == "sparse" && region_mask(rand_z, rand_x) == 1)
                    
                    indices_to_remove = [indices_to_remove; j + prev_size]; % Mark this index for removal

                end
            end
        end
        total_num_fn = round(length(indices_to_remove)* fn_rate);
        % Randomly remove points from this frame if they are in the desired region
        if ~isempty(indices_to_remove)
            num_to_remove = min(length(indices_to_remove), total_num_fn);
            remove_idx = indices_to_remove(randperm(length(indices_to_remove), num_to_remove));
            
            % Remove the selected points from the current frame
            data_with_fn(remove_idx, :) = [];
        else
            display("not enough points in the region")
        end
        
    end
end

function data_with_fn = addFalseNegativesTotalNumberSpecificRegion(data, num_fn, region_mask, region_type, x_grid, z_grid)
    % Adds a specific number of false negatives to either dense or sparse regions.
    % num_fn: Total number of false negatives to introduce.
    % region_mask: Binary mask where 1 indicates a region of interest (dense/sparse), 0 otherwise.
    % region_type: Either 'dense' or 'sparse'

    data_with_fn = data;
    frames = unique(data(:, 1));  % Unique frame numbers
    total_points = size(data, 1);
    
    % Ensure we don't try to remove more points than we have
    num_fn = min(num_fn, total_points);
    
    points_removed = 0; % Track how many points have been removed
    indices_to_remove = [];
    prev_size = 0;
    if num_fn ~= 0
        for i = 1:length(frames)
            if points_removed >= num_fn
                break;
            end
            
            frame_idx = frames(i);
            frame_data = data(data(:, 1) == frame_idx, :);
            
            % Extract the true positive coordinates (x, z) for the current frame
            true_x = frame_data(:, 3); % Assuming x-coordinates are in the 3rd column
            true_z = frame_data(:, 5); % Assuming z-coordinates are in the 5th column
            
            % Check if the points are in the desired region
            prev_size = prev_size + size(data(data(:, 1) == frame_idx-1, :),1);
            for j = 1:size(frame_data, 1)
                % Convert true positive (x, z) to grid coordinates
                [~, rand_x] = min(abs(x_grid(1, :) - true_x(j)));
                [~, rand_z] = min(abs(z_grid(:, 1) - true_z(j)));
                
                % Check if the point is in the desired region (dense or sparse)
                if (region_type == "dense" && region_mask(rand_z, rand_x) == 1) || ...
                   (region_type == "sparse" && region_mask(rand_z, rand_x) == 1)
                    indices_to_remove = [indices_to_remove; j + prev_size]; % Mark this index for removal
                end
            end
        end
        
        % Remove the desired number of points
        if ~isempty(indices_to_remove)
            num_to_remove = min(length(indices_to_remove), num_fn);
            remove_idx = indices_to_remove(randperm(length(indices_to_remove), num_to_remove));
            
            % Remove the selected points from the current frame
            data_with_fn(remove_idx, :) = [];
        else
            display("not enough points in the region")
        end
        
    else
        data_with_fn = data;
    end
end

function data_with_fn = addFalseNegativesPerFrameNumberSpecificRegion(data, num_fn_per_frame, region_mask, region_type, x_grid, z_grid)
    % Adds a specific number of false negatives per frame to either dense or sparse regions.
    % num_fn_per_frame: Number of false negatives to introduce per frame.
    % region_mask: Binary mask where 1 indicates a region of interest (dense/sparse), 0 otherwise.
    % region_type: Either 'dense' or 'sparse'

    data_with_fn = [];
    frames = unique(data(:, 1));  % Unique frame numbers
    
    if num_fn_per_frame ~= 0
        for i = 1:length(frames)
            frame_idx = frames(i);
            frame_data = data(data(:, 1) == frame_idx, :);
            
            % Extract the true positive coordinates (x, z) for the current frame
            true_x = frame_data(:, 3); % Assuming x-coordinates are in the 3rd column
            true_z = frame_data(:, 5); % Assuming z-coordinates are in the 5th column
            
            % Find points that lie within the desired region
            indices_to_remove = [];
            for j = 1:size(frame_data, 1)
                % Convert true positive (x, z) to grid coordinates
                [~, rand_x] = min(abs(x_grid(1, :) - true_x(j)));
                [~, rand_z] = min(abs(z_grid(:, 1) - true_z(j)));
                
                % Check if the point is in the desired region (dense or sparse)
                if (region_type == "dense" && region_mask(rand_z, rand_x) == 1) || ...
                   (region_type == "sparse" && region_mask(rand_z, rand_x) == 1)
                    indices_to_remove = [indices_to_remove; j ]; % Mark this index for removal
                end
            end
            if ~isempty(indices_to_remove)
                num_to_remove = min(length(indices_to_remove), num_fn_per_frame);
                remove_idx = indices_to_remove(randperm(length(indices_to_remove), num_to_remove));
                frame_data(remove_idx,:) = [];
                data_with_fn = [data_with_fn; frame_data];
            else
                disp("not enough data in the region to remove")
            end
            
        end
    else
        data_with_fn = data;
    end
end

function modified_data = modifyDataPerFrameNumberSpecificRegion(data, x_grid, z_grid, num_fn_per_frame, num_fp_per_frame, threshold, region_mask, region_type)
    % Step 1: Add false negatives to each frame based on the region type
    data_with_fn = addFalseNegativesPerFrameNumberSpecificRegion(data, num_fn_per_frame, region_mask, region_type, x_grid, z_grid);
    
    % Step 2: Add false positives to each frame based on the region type
    modified_data = addFalsePositivesPerFrameNumberSpecificRegion(data_with_fn, x_grid, z_grid, num_fp_per_frame, threshold, region_mask, region_type);
end


function modified_data = modifyDataTotalNumberSpecificRegion(data, x_grid, z_grid, num_fn_total, num_fp_total, threshold, region_mask, region_type)
    % Step 1: Add false negatives to all frames based on the region type
    data_with_fn = addFalseNegativesTotalNumberSpecificRegion(data, num_fn_total, region_mask, region_type, x_grid, z_grid);
    
    % Step 2: Add false positives to all frames based on the region type
    modified_data = addFalsePositivesTotalNumberSpecificRegion(data_with_fn, x_grid, z_grid, num_fp_total, threshold, region_mask, region_type);
end


function modified_data = modifyDataTotalRateSpecificRegion(data, x_grid, z_grid, fn_rate, fp_rate, threshold, region_mask, region_type)
    % Step 1: Add false negatives to the total data based on the region type
    data_with_fn = addFalseNegativesTotalRateSpecificRegion(data, fn_rate, region_mask, region_type, x_grid, z_grid);
    
    % Step 2: Calculate the number of false positives to add
    total_points = size(data, 1);
    num_to_add = round(total_points * fp_rate);
    
    % Step 3: Add false positives to the data with false negatives based on the region type
    modified_data = addFalsePositivesTotalNumberSpecificRegion(data_with_fn, x_grid, z_grid, num_to_add, threshold, region_mask, region_type);
end


function modified_data = modifyDataPerFrameRateSpecificRegion(data, x_grid, z_grid, fn_rate, fp_rate, threshold, region_mask, region_type)
    % Step 1: Add false negatives to each frame based on the region type
    data_with_fn = addFalseNegativesPerFrameRateSpecificRegion(data, fn_rate, region_mask, region_type, x_grid, z_grid);
    
    % Step 2: Calculate the number of false positives to add
    total_points = size(data, 1);
    num_to_add = round(total_points * fp_rate);
    
    % Step 3: Add false positives to each frame based on the region type
    modified_data = addFalsePositivesPerFrameRateSpecificRegion(data_with_fn, data, x_grid, z_grid, fp_rate, threshold, region_mask, region_type);
end

function plot_SR_map(SR_density_resized,titleC)

    figure;
    imagesc([-0.015 0.015], [0.047 0.107], SR_density_resized', [0 5]);
    colormap('hot');
    axis image;
    % Remove the numbers on the x and y-axis
    % xticks([]);
    % yticks([]);
    % Add a white line on the left below part of the image with a desired length (e.g., 2 mm)
    line([-0.014, -0.013], [0.105, 0.105], 'Color', 'white', 'LineWidth', 1);
    % Add the text "1mm" above the line in white
    text(-0.015, 0.103, '1mm', 'Color', 'white', 'FontSize', 12);
    title(titleC)
    xticks([]);
    yticks([]);
end


function ssim_value = calculate_ssim(SR_density_gt_resized, SR_density_mod_resized, PxSet)
    % Compute the SR density maps

    % Compute SSIM
    ssim_value = ssim(SR_density_mod_resized, SR_density_gt_resized);
end

function dice_coeff = calculate_dice(SR_density_gt_resized, SR_density_mod_resized, PxSet, threshold)

    % Threshold the density maps to create binary maps
    binary_gt = SR_density_gt_resized > threshold;
    binary_mod = SR_density_mod_resized > threshold;

    % Compute Dice coefficient
    intersection = sum(binary_gt(:) & binary_mod(:));
    union = sum(binary_gt(:)) + sum(binary_mod(:));
    
    dice_coeff = 2 * intersection / union;
end

function psnr_value = calculate_psnr(SR_density_gt_resized, SR_density_mod_resized, PxSet)

    % Compute Mean Squared Error (MSE)
    mse_value = mean((SR_density_gt_resized(:) - SR_density_mod_resized(:)).^2);

    % Compute PSNR
    max_pixel_value = max(SR_density_gt_resized(:));
    psnr_value = 10 * log10(max_pixel_value^2 / mse_value);
end


function region_type_map = classifyRegions(kde_density, threshold, region_type)
    if strcmp(region_type, 'dense')
        % Dense regions are those where KDE density is greater than the threshold
        region_type_map = kde_density > threshold;
    else
        % Sparse regions are those where KDE density is less than the threshold
        region_type_map = kde_density <= threshold;
    end
end

function ssim_value = calculate_ssim_mask(SR_density_gt_resized, SR_density_mod_resized, mask)
    % Apply the mask to focus on specific regions
    SR_density_gt_masked = SR_density_gt_resized(mask);
    SR_density_mod_masked = SR_density_mod_resized(mask);

    % Compute SSIM for the masked regions
    ssim_value = ssim(SR_density_mod_masked, SR_density_gt_masked);
end

function dice_coeff = calculate_dice_mask(SR_density_gt_resized, SR_density_mod_resized, mask)
    % Threshold and apply the mask to create binary maps
    binary_gt = (SR_density_gt_resized > 0.0) & mask;
    binary_mod = (SR_density_mod_resized > 0.0) & mask;

    % Compute Dice coefficient for the masked regions
    intersection = sum(binary_gt(:) & binary_mod(:));
    union = sum(binary_gt(:)) + sum(binary_mod(:));
    
    dice_coeff = 2 * intersection / union;
end

function psnr_value = calculate_psnr_mask(SR_density_gt_resized, SR_density_mod_resized, mask)
    % Apply the mask to focus on specific regions
    SR_density_gt_masked = SR_density_gt_resized(mask);
    SR_density_mod_masked = SR_density_mod_resized(mask);

    % Compute Mean Squared Error (MSE) for the masked regions
    mse_value = mean((SR_density_gt_masked(:) - SR_density_mod_masked(:)).^2);

    % Compute PSNR for the masked regions
    max_pixel_value = max(SR_density_gt_masked(:));
    psnr_value = 10 * log10(max_pixel_value^2 / mse_value);
end
