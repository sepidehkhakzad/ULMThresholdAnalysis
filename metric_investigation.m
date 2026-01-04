
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Ultrasound Super-Resolution Imaging: Threshold Analysis Tool
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
% quality metrics (SSIM, Dice coefficient, PSNR) under various error
% scenarios.
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

% Load the ground truth data
ground_truth = load('gt_scat_inside2.txt');

% Load the metadata file
load('metadata_simu2.mat');
x_grid = repmat((PxSet.mapX), 512, 1);
z_grid = repmat((PxSet.mapZ)', 1, 512);

cf = SimSet.centre_frequency;
% Calculate the metrics
lambdaa = 1540/cf;
threshold = lambdaa/2;   % Example threshold for matching
gt_SRDensity = compute_gt_density(ground_truth, PxSet);


% Define the ranges for num_fp_per_frame and fn_rate
num_fp_range = 0:10:100;  % range for number of false positives per frame
fn_rate_range = 0:0.1:0.8;  % range for false negative rate
fp_rate_range = 0:0.1:0.8;
num_fn_range = 0:10:100;

% Initialize matrices to store the metrics
ssim_matrix_number = zeros(length(num_fp_range), length(num_fn_range), 2);
dice_matrix_number = zeros(length(num_fp_range), length(num_fn_range), 2);
psnr_matrix_number = zeros(length(num_fp_range), length(num_fn_range), 2);

ssim_matrix_rate = zeros(length(fp_rate_range), length(fn_rate_range), 2);
dice_matrix_rate = zeros(length(fp_rate_range), length(fn_rate_range), 2);
psnr_matrix_rate = zeros(length(fp_rate_range), length(fn_rate_range), 2);

ssim_tmp = cell(1,4);
dice_tmp = cell(1,4);
psnr_tmp = cell(1,4);

% mode = 1: Total Number
% mode = 2: Per-Frame Number
% mode = 3: Total Rate
% mode = 4: Per-Frame Rate
% This will take multiple hours to run.
parfor mode = 1:4

    isNumber   = mode <= 2;
    isPerFrame = mod(mode,2) == 0;

    if isNumber
        fp_range = num_fp_range;
        fn_range = num_fn_range;

        ssim_local = zeros(length(fp_range), length(fn_range));
        dice_local = zeros(length(fp_range), length(fn_range));
        psnr_local = zeros(length(fp_range), length(fn_range));
    else
        fp_range = fp_rate_range;
        fn_range = fn_rate_range;

        ssim_local = zeros(length(fp_range), length(fn_range));
        dice_local = zeros(length(fp_range), length(fn_range));
        psnr_local = zeros(length(fp_range), length(fn_range));
    end

    for i = 1:length(fp_range)
        for j = 1:length(fn_range)

            fp_val = fp_range(i);
            fn_val = fn_range(j);

            if isNumber && ~isPerFrame
                modified_data = modifyDataTotalNumber( ...
                    ground_truth, x_grid, z_grid, fp_val, fn_val, threshold);

            elseif isNumber && isPerFrame
                modified_data = modifyDataPerFrameNumber( ...
                    ground_truth, x_grid, z_grid, fp_val, fn_val, threshold);

            elseif ~isNumber && ~isPerFrame
                modified_data = modifyDataTotalRate( ...
                    ground_truth, x_grid, z_grid, fp_val, fn_val, threshold);

            else % ~isNumber && isPerFrame
                modified_data = modifyDataPerFrameRate( ...
                    ground_truth, x_grid, z_grid, fp_val, fn_val, threshold);
            end

            modified_SRDensity = compute_gt_density(modified_data, PxSet);

            ssim_local(i,j) = calculate_ssim(gt_SRDensity, modified_SRDensity, PxSet);
            dice_local(i,j) = calculate_dice(gt_SRDensity, modified_SRDensity, PxSet);
            psnr_local(i,j) = calculate_psnr(gt_SRDensity, modified_SRDensity, PxSet);
        end
    end

    ssim_tmp{mode} = ssim_local;
    dice_tmp{mode} = dice_local;
    psnr_tmp{mode} = psnr_local;
end

% Number-based
ssim_matrix_number(:,:,1) = ssim_tmp{1};
ssim_matrix_number(:,:,2) = ssim_tmp{2};
dice_matrix_number(:,:,1) = dice_tmp{1};
dice_matrix_number(:,:,2) = dice_tmp{2};
psnr_matrix_number(:,:,1) = psnr_tmp{1};
psnr_matrix_number(:,:,2) = psnr_tmp{2};

% Rate-based
ssim_matrix_rate(:,:,1) = ssim_tmp{3};
ssim_matrix_rate(:,:,2) = ssim_tmp{4};
dice_matrix_rate(:,:,1) = dice_tmp{3};
dice_matrix_rate(:,:,2) = dice_tmp{4};
psnr_matrix_rate(:,:,1) = psnr_tmp{3};
psnr_matrix_rate(:,:,2) = psnr_tmp{4};


% Create meshgrid for plotting
[fp_rate_range_grid, fn_rate_grid] = meshgrid(fp_rate_range, fn_rate_range);
% Plot SSIM
figure;
surf(fp_rate_range_grid, fn_rate_grid, ssim_matrix_rate(:,:,1));
xlabel('False Negative Rate Total');
ylabel('False Positive Rate Total');
zlabel('SSIM');
title('SSIM vs. False Negative Rate and False Positive Rate Total');
colorbar;
figure;
surf(fp_rate_range_grid, fn_rate_grid, dice_matrix_rate(:,:,1));
xlabel('False Negative Rate Total');
ylabel('False Positive Rate Total');
zlabel('Dice');
title('Dice vs. False Negative Rate and False Positive Rate Total');
colorbar;
figure;
surf(fp_rate_range_grid, fn_rate_grid, psnr_matrix_rate(:,:,1));
xlabel('False Negative Rate Total');
ylabel('False Positive Rate Total');
zlabel('PSNR');
title('PSNR vs. False Negative Rate and False Positive Rate Total');
colorbar;
% Plot SSIM
figure;
surf(fp_rate_range_grid, fn_rate_grid, ssim_matrix_rate(:,:,2));
xlabel('False Negative Rate Per Frame');
ylabel('False Positive Rate Per Frame');
zlabel('SSIM');
title('SSIM vs. False Negative Rate and False Positive Rate Per Frame');
colorbar;
figure;
surf(fp_rate_range_grid, fn_rate_grid, dice_matrix_rate(:,:,2));
xlabel('False Negative Rate Per Frame');
ylabel('False Positive Rate Per Frame');
zlabel('Dice');
title('Dice vs. False Negative Rate and False Positive Rate Per Frame');
colorbar;
figure;
surf(fp_rate_range_grid, fn_rate_grid, psnr_matrix_rate(:,:,2));
xlabel('False Negative Rate Per Frame');
ylabel('False Positive Rate Per Frame');
zlabel('PSNR');
title('PSNR vs. False Negative Rate and False Positive Rate Per Frame');
colorbar;
% Plot SSIM
[num_fp_range_grid, num_fn_range_grid] = meshgrid(num_fp_range, num_fn_range);
figure;
surf(num_fp_range_grid, num_fn_range_grid, ssim_matrix_number(:,:,2));
xlabel('False Negative Numbers Per Frame');
ylabel('False Positive Numbers Per Frame');
zlabel('SSIM');
title('SSIM vs. False Negative Numbers and False Positive Numbers Per Frame');
colorbar;
figure;
surf(num_fp_range_grid, num_fn_range_grid, dice_matrix_number(:,:,2));
xlabel('False Negative Numbers Per Frame');
ylabel('False Positive Numbers Per Frame');
zlabel('Dice');
title('Dice vs. False Negative Numbers and False Positive Numbers Per Frame');
colorbar;
figure;
surf(num_fp_range_grid, num_fn_range_grid, psnr_matrix_number(:,:,2));
xlabel('False Negative Numbers Per Frame');
ylabel('False Positive Numbers Per Frame');
zlabel('PSNR');
title('PSNR vs. False Negative Numbers and False Positive Numbers Per Frame');
colorbar;

% Plot SSIM
[num_fp_range_grid, num_fn_range_grid] = meshgrid(num_fp_range, num_fn_range);
figure;
surf(num_fp_range_grid, num_fn_range_grid, ssim_matrix_number(:,:,2));
xlabel('False Negative Numbers Per Frame');
ylabel('False Positive Numbers Per Frame');
zlabel('SSIM');
title('SSIM vs. False Negative Numbers and False Positive Numbers Per Frame');
colorbar;
figure;
surf(num_fp_range_grid, num_fn_range_grid, dice_matrix_number(:,:,2));
xlabel('False Negative Numbers Per Frame');
ylabel('False Positive Numbers Per Frame');
zlabel('Dice');
title('Dice vs. False Negative Numbers and False Positive Numbers Per Frame');
colorbar;
figure;
surf(num_fp_range_grid, num_fn_range_grid, psnr_matrix_number(:,:,2));
xlabel('False Negative Numbers Per Frame');
ylabel('False Positive Numbers Per Frame');
zlabel('PSNR');
title('PSNR vs. False Negative Numbers and False Positive Numbers Per Frame');
colorbar;
figure;
surf(num_fp_range_grid, num_fn_range_grid, ssim_matrix_number(:,:,1));
xlabel('False Negative Numbers Total');
ylabel('False Positive Numbers Total');
zlabel('SSIM');
title('SSIM vs. False Negative Numbers and False Positive Numbers Total');
colorbar;
figure;
surf(num_fp_range_grid, num_fn_range_grid, dice_matrix_number(:,:,1));
xlabel('False Negative Numbers Total');
ylabel('False Positive Numbers Total');
zlabel('Dice');
title('Dice vs. False Negative Numbers and False Positive Numbers Total');
colorbar;
figure;
surf(num_fp_range_grid, num_fn_range_grid, psnr_matrix_number(:,:,1));
xlabel('False Negative Numbers Total');
ylabel('False Positive Numbers Total');
zlabel('PSNR');
title('PSNR vs. False Negative Numbers and False Positive Numbers Total');
colorbar;


function data_with_fp = addFalsePositivesPerFrameRate(data, x_grid, z_grid, fp_percentage, distance_threshold)
    data_with_fp = data;
    frames = unique(data(:, 1)); % Unique frame numbers
    
    for i = 1:length(frames)
        frame_idx = frames(i);
        frame_data = data(data(:, 1) == frame_idx, :);
        
        % Extract the true positive coordinates (x, z) for the current frame
        true_x = frame_data(:, 3); % Assuming x-coordinates are in the 3rd column
        true_z = frame_data(:, 5); % Assuming z-coordinates are in the 5th column
        
        % Calculate the number of false positives to add for this frame
        num_frame_points = size(frame_data, 1);
        num_fp_for_frame = round(num_frame_points * fp_percentage / 100);
        
        for j = 1:num_fp_for_frame
            false_positive_added = false;
            while ~false_positive_added
                % Choose random position in the grid
                rand_x = randi(size(x_grid, 2)); % Random column index
                rand_z = randi(size(z_grid, 1)); % Random row index
                
                % Get the x and z coordinates of the randomly selected grid point
                fp_x = x_grid(rand_z, rand_x);
                fp_z = z_grid(rand_z, rand_x);
                
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


function data_with_fp = addFalsePositivesTotalRate(data, x_grid, z_grid, fp_percentage, distance_threshold)
    % Calculate the number of false positives to add
    num_original_points = size(data, 1);
    total_num_fp = round(num_original_points * fp_percentage / 100);
    
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


function data_with_fp = addFalsePositivesTotalNumber(data, x_grid, z_grid, total_num_fp, distance_threshold)
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

function data_with_fp = addFalsePositivesPerFrameNumber(data, x_grid, z_grid, num_fp_per_frame, distance_threshold)
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
                rand_x = randi(size(x_grid, 2));  % Random column index
                rand_z = randi(size(z_grid, 1));  % Random row index
                
                % Get the x and z coordinates of the randomly selected grid point
                fp_x = x_grid(rand_z, rand_x);
                fp_z = z_grid(rand_z, rand_x);
                
                % Compute the distance to all true positives in this frame
                distances = sqrt((true_x - fp_x).^2 + (true_z - fp_z).^2);
                
                % If the minimum distance is greater than the threshold, accept this point
                if all(distances > distance_threshold)
                    % Create a random false positive
                    false_positive = [frame_idx, max(frame_data(:, 2)) + 1, fp_x, 0, fp_z, 0];
                    
                    % Add the false positive to the data
                    data_with_fp = [data_with_fp; false_positive];
                    false_positive_added = true;
                end
            end
        end
    end
end


function data_with_fn = addFalseNegativesPerFrameRate(data, fn_rate)
    data_with_fn = [];
    frames = unique(data(:, 1));  % Unique frame numbers
    
    for i = 1:length(frames)
        frame_idx = frames(i);
        frame_data = data(data(:, 1) == frame_idx, :);
        
        % Calculate the number of false negatives to introduce
        num_to_remove = round(size(frame_data, 1) * fn_rate);
        
        % Randomly select rows to remove
        remove_idx = randperm(size(frame_data, 1), num_to_remove);
        frame_data(remove_idx, :) = [];
        
        data_with_fn = [data_with_fn; frame_data];
    end
end

function data_with_fn = addFalseNegativesTotalRate(data, fn_rate)
    % Calculate the total number of points to remove
    total_points = size(data, 1);
    num_to_remove = round(total_points * fn_rate);
    
    % Randomly select indices to remove
    all_indices = 1:total_points;
    remove_indices = randsample(all_indices, num_to_remove, false);
    
    % Remove the selected points
    data_with_fn = data;
    data_with_fn(remove_indices, :) = [];
end


function data_with_fn = addFalseNegativesTotalNumber(data, num_fn)
    total_points = size(data, 1);
    
    % Ensure we're not trying to remove more points than exist
    num_fn = min(num_fn, total_points);
    
    % Randomly select indices to remove
    all_indices = 1:total_points;
    remove_indices = randsample(all_indices, num_fn, false);
    
    % Remove the selected points
    data_with_fn = data;
    data_with_fn(remove_indices, :) = [];
end

function data_with_fn = addFalseNegativesPerFrameNumber(data, num_fn_per_frame)
    data_with_fn = [];
    frames = unique(data(:, 1)); % Unique frame numbers
    
    for i = 1:length(frames)
        frame_idx = frames(i);
        frame_data = data(data(:, 1) == frame_idx, :);
        
        % Ensure we're not trying to remove more points than exist in this frame
        num_fn = min(num_fn_per_frame, size(frame_data, 1));
        
        if num_fn > 0
            % Randomly select rows to remove
            remove_idx = randsample(size(frame_data, 1), num_fn, false);
            frame_data(remove_idx, :) = [];
        end
        
        data_with_fn = [data_with_fn; frame_data];
    end
end

function [metrics] = calculateMetrics(ground_truth, modified_data, x_grid, z_grid, threshold)
    metrics = struct('TP', 0, 'FP', 0, 'FN', 0, 'TN', 0, 'precision', 0, 'recall', 0, 'specificity', 0, 'sensitivity', 0);
    
    frames = unique(ground_truth(:, 1));  % Unique frame numbers
    
    for i = 1:length(frames)
        frame_idx = frames(i);
        
        % Extract data for the current frame
        gt_frame_data = ground_truth(ground_truth(:, 1) == frame_idx, :);
        mod_frame_data = modified_data(modified_data(:, 1) == frame_idx, :);
        
        % Initialize counters for this frame
        TP = 0; FP = 0; FN = 0; TN = 0;
        
        % Compare ground truth with modified data
        for j = 1:size(gt_frame_data, 1)
            gt_x = gt_frame_data(j, 3);
            gt_z = gt_frame_data(j, 5);
            
            % Find the closest match in the modified data
            distances = sqrt((mod_frame_data(:, 3) - gt_x).^2 + (mod_frame_data(:, 5) - gt_z).^2);
            [min_distance, min_idx] = min(distances);
            
            if min_distance <= threshold
                TP = TP + 1;  % True Positive
                mod_frame_data(min_idx, :) = [];  % Remove matched item to avoid multiple counting
            else
                FN = FN + 1;  % False Negative
            end
        end
        
        % Remaining items in modified data are False Positives
        FP = FP + size(mod_frame_data, 1);
        
        % TN is not easily definable in this context, so it's left as 0
        
        % Aggregate the frame results
        metrics.TP = metrics.TP + TP;
        metrics.FP = metrics.FP + FP;
        metrics.FN = metrics.FN + FN;
        % TN would normally require a true negative area calculation, which is complex without additional context
    end
    
    % Calculating Precision, Recall, Specificity, and Sensitivity
    metrics.precision = metrics.TP / (metrics.TP + metrics.FP);
    metrics.recall = metrics.TP / (metrics.TP + metrics.FN);
    % Specificity and Sensitivity require TN, which is not computed here. 
    % We assume specificity = recall, and sensitivity is based on TP / (TP + FN)
    metrics.specificity = metrics.recall;  % Placeholder without TN calculation
    metrics.sensitivity = metrics.recall;
end

function modified_data = modifyDataPerFrameNumber(data, x_grid, z_grid, num_fn_per_frame, num_fp_per_frame, threshold)
    % Step 1: Add false negatives to each frame
    data_with_fn = addFalseNegativesPerFrameNumber(data, num_fn_per_frame);
    
    % Step 2: Add false positives to each frame
    modified_data = addFalsePositivesPerFrameNumber(data_with_fn, x_grid, z_grid, num_fp_per_frame, threshold);
end

function modified_data = modifyDataTotalNumber(data, x_grid, z_grid, num_fn_per_frame, num_fp_per_frame, threshold)
    % Step 1: Add false negatives to all frames
    data_with_fn = addFalseNegativesTotalNumber(data, num_fn_per_frame);
    
    % Step 2: Add false positives to all frames
    modified_data = addFalsePositivesTotalNumber(data_with_fn, x_grid, z_grid, num_fp_per_frame,threshold);
end

function modified_data = modifyDataTotalRate(data, x_grid, z_grid, fn_rate, fp_rate, threshold)
    % Step 1: Add false negatives to the total data
    data_with_fn = addFalseNegativesTotalRate(data, fn_rate);
    total_points = size(data, 1);
    num_to_add = round(total_points * fp_rate);
    % Step 3: Add false positives to the data with false negatives
    modified_data = addFalsePositivesTotalNumber(data_with_fn, x_grid, z_grid, num_to_add,threshold);
end

function modified_data = modifyDataPerFrameRate(data, x_grid, z_grid, fn_rate, fp_rate, threshold)
    % Step 1: Add false negatives to each frame
    data_with_fn = addFalseNegativesPerFrameRate(data, fn_rate);
    total_points = size(data, 1);
    num_to_add = round(total_points * fp_rate);
    % Step 2: Add false positives to each frame
    modified_data = addFalsePositivesPerFrameNumber(data_with_fn, x_grid, z_grid, num_to_add, threshold);
end


function SR_density_resized = compute_gt_density(gt, PxSet)
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

function dice_coeff = calculate_dice(SR_density_gt_resized, SR_density_mod_resized, PxSet)

    % Threshold the density maps to create binary maps
    binary_gt = SR_density_gt_resized > 0.0;
    binary_mod = SR_density_mod_resized > 0.0;

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



